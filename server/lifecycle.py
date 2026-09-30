"""Quit completely when the app's window is closed.

Clone Bench is a local server with a browser tab as its window. Left alone the server
keeps running after the tab is closed, which is how a machine ends up with several idle
copies listening on the port. This module makes the window and the server end together:

* An open page holds a connection to the server (the window channel) and also sends a
  heartbeat every few seconds. Either one is enough to say somebody is looking at the app.
* When that connection drops, or the page says it is closing, the window has gone away, and
  the server stops six seconds later unless a page connects again. A reload reconnects in
  about a second, so reloading does not quit the app. The dropped connection is the signal
  that matters: a browser is free to close a tab without telling the page, and some do.
* If no page has been connected for three minutes and nothing said it was closing (the
  browser crashed, or was force-quit), the server stops as well.
* If the launcher that started the server is gone — the app was quit from the Dock, which
  is a SIGKILL a shell trap never sees — the server stops.
* Whatever the reason, a process the app started is killed with it, as a group.

Counting is in checks rather than wall-clock time, so a Mac waking from sleep does not
find the deadline already past and quit the app the user just came back to.
"""
from __future__ import annotations

import atexit
import os
import signal
import sys
import threading
import time

CHECK_S = 2.0            # how often the server looks at who is connected
STALE_S = 8.0            # a page that has not sent a heartbeat for this long is gone (it sends every 3 s)
AFTER_CLOSE = 3          # checks with nobody connected after a page said it was closing (about 6 s)
AFTER_SILENCE = 90       # checks with nobody connected and no such message (about 3 min)
LAUNCHER_ENV = "CB_LAUNCHER_PID"
DISABLE_ENV = "CB_NO_AUTOQUIT"    # set by dev.sh: a development server stays up on its own

_lock = threading.Lock()
_pages: dict[str, float] = {}       # page id -> when it last said it was there
_windows: dict[int, str] = {}       # open window channel -> the page holding it (see /api/window)
_state = {"started": False, "closing": False, "children": [], "next": 0}


def _log(text: str) -> None:
    print(f"[{time.strftime('%F %T')}] lifecycle: {text}", file=sys.stderr, flush=True)


def _alive(pid: int) -> bool:
    try:
        os.kill(pid, 0)
    except ProcessLookupError:
        return False
    except PermissionError:
        return True
    return True


# ---------------------------------------------------------------- what the page reports
def ping(page_id: str) -> None:
    """A page is open and looking at the app."""
    with _lock:
        _pages[page_id] = time.monotonic()


def closing(page_id: str) -> None:
    """A page is going away. If it was the last one, the app has been closed."""
    with _lock:
        _pages.pop(page_id, None)
        _state["closing"] = True


def window_opened(page_id: str = "") -> int:
    """A page has opened its window channel. Returns the token to close it with."""
    with _lock:
        _state["next"] += 1
        token = _state["next"]
        _windows[token] = page_id
        return token


def window_closed(token: int) -> None:
    """The channel dropped: the window is gone, whether or not the page got to say so.

    The page's heartbeat goes with it. Leaving it to go stale on its own would keep the app
    looking occupied for several more seconds, and every one of those seconds would clear the
    flag this just set — the app would then wait the full three minutes instead of six seconds.
    """
    with _lock:
        page_id = _windows.pop(token, "")
        _pages.pop(page_id, None)
        _state["closing"] = True


def connected() -> int:
    """How many pages are looking at the app. A page holds a window channel and sends a
    heartbeat, so it is counted once, and either signal alone is enough to keep it alive."""
    now = time.monotonic()
    with _lock:
        for pid in [k for k, t in _pages.items() if now - t > STALE_S]:
            _pages.pop(pid, None)
        return max(len(_windows), len(_pages))


# ---------------------------------------------------------------- processes the app started
def register_child(proc) -> None:
    """A process started by the app (in its own process group) that must not outlive it."""
    with _lock:
        _state["children"] = [p for p in _state["children"] if p.poll() is None] + [proc]


def stop_children() -> None:
    with _lock:
        children, _state["children"] = _state["children"], []
    for proc in children:
        if proc.poll() is None:
            try:
                os.killpg(proc.pid, signal.SIGTERM)   # the child and anything it started share one group
            except (ProcessLookupError, PermissionError, OSError):
                pass


# ---------------------------------------------------------------- the rule
def should_quit(connected: int, seen: bool, empty: int, closing: bool) -> tuple[bool, bool, int, bool]:
    """One step of the rule. Returns (quit now?, seen, empty checks, closing)."""
    if connected > 0:
        return False, True, 0, False
    if not seen:                                      # nobody has opened the app yet
        return False, False, 0, closing
    empty += 1
    return (closing and empty >= AFTER_CLOSE) or empty >= AFTER_SILENCE, True, empty, closing


def quit_now(reason: str) -> None:
    _log(f"quitting: {reason}")
    stop_children()
    os._exit(0)


def _watch() -> None:
    launcher = int(os.environ.get(LAUNCHER_ENV, "0") or 0)
    seen, empty = False, 0
    while True:
        time.sleep(CHECK_S)
        if launcher and not _alive(launcher):
            quit_now("the app was quit")
        quit_it, seen, empty, closes = should_quit(connected(), seen, empty, _state["closing"])
        _state["closing"] = closes
        if quit_it:
            quit_now("its window was closed" if closes else "no window has been connected for three minutes")


def start() -> None:
    """Start watching (once per server)."""
    with _lock:
        if _state["started"]:
            return
        _state["started"] = True
    atexit.register(stop_children)
    if os.environ.get(DISABLE_ENV):
        _log("automatic quit is off (development server)")
        return
    threading.Thread(target=_watch, daemon=True).start()
