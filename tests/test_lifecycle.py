"""The rule that decides when Clone Bench quits, tested without starting a server.

Run:  python3 tests/test_lifecycle.py
"""
import os
import sys
import time

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from server import lifecycle as L  # noqa: E402

fails = []


def check(name, got, want):
    if got != want:
        fails.append(f"{name}: got {got!r}, wanted {want!r}")
    print(f"{'ok  ' if got == want else 'FAIL'} {name}")


def run(steps, seen=False, empty=0, closing=False):
    """Feed a list of 'how many pages are connected' through the rule; return the check it quits on."""
    for i, connected in enumerate(steps, 1):
        quit_it, seen, empty, closing = L.should_quit(connected, seen, empty, closing)
        if quit_it:
            return i
    return None


# nobody has ever opened the app: it waits, however long that takes
check("an app nobody has opened never quits", run([0] * 500), None)

# a page is open, then the window is closed: the beacon sets closing, and it quits about 6 s later
check("window closed -> quits on the third empty check",
      run([0] * 10, seen=True, closing=True), L.AFTER_CLOSE)

# a reload: the page goes away, says so, and is back before the third check
check("a reload does not quit it", run([1, 0, 0, 1, 1, 1], closing=True), None)

# after that reload the app must be willing to quit again when the new window is closed
seen, empty, closing = True, 0, True
for connected in [1, 0, 0, 1]:                      # reload
    _, seen, empty, closing = L.should_quit(connected, seen, empty, closing)
check("the closing flag is cleared by a page that connects", closing, False)

# browser crashed or force-quit: no beacon, so it waits three minutes
check("a silent disconnect waits for AFTER_SILENCE", run([0] * 200, seen=True), L.AFTER_SILENCE)
check("three minutes is about right", round(L.AFTER_SILENCE * L.CHECK_S), 180)
check("six seconds is about right", round(L.AFTER_CLOSE * L.CHECK_S), 6)

# the heartbeat bookkeeping
L._pages.clear()
L.ping("a")
L.ping("b")
check("two open pages count as two", L.connected(), 2)
L.closing("a")
check("a page that said it is closing is gone", L.connected(), 1)
L._pages["b"] = time.monotonic() - (L.STALE_S + 1)
check("a page that stopped its heartbeat is gone", L.connected(), 0)
check("the closing flag was set by the beacon", L._state["closing"], True)

# a heartbeat interval that does not fit inside the staleness window would drop live pages
check("the page's 3 s heartbeat fits inside the staleness window", L.STALE_S > 3 * 2, True)

# The window channel is the signal that a window has gone; the heartbeat is only a fallback. When a
# channel drops, its page's heartbeat has to go with it. Left to go stale on its own it keeps the app
# looking occupied, and the watcher clears the "closing" flag on every check where somebody is
# connected — which is how a closed window turned into a three-minute wait instead of six seconds.
L._pages.clear(); L._windows.clear(); L._state["closing"] = False
tok = L.window_opened("p1")
L.ping("p1")
check("a page with a channel counts once", L.connected(), 1)
L.window_closed(tok)
check("closing the channel takes its heartbeat with it", L.connected(), 0)
check("and marks the app as closed", L._state["closing"], True)

# the watcher clears that flag while anyone is connected; with the fix there is nobody left to clear it
seen, empty, closing = True, 0, L._state["closing"]
quit_it, seen, empty, closing = L.should_quit(L.connected(), seen, empty, closing)
check("the flag survives the next check", closing, True)

# two windows: closing one does not close the app
L._pages.clear(); L._windows.clear(); L._state["closing"] = False
a, b = L.window_opened("a"), L.window_opened("b")
L.ping("a"); L.ping("b")
L.window_closed(a)
check("a second window keeps the app open", L.connected(), 1)
L.window_closed(b)
check("closing the last window closes the app", L.connected(), 0)

# a process the app started must not outlive it, and nor must anything it started
import subprocess  # noqa: E402

parent = subprocess.Popen([sys.executable, "-c",
                           "import subprocess,sys,time;"
                           "subprocess.Popen([sys.executable,'-c','import time;time.sleep(60)']);"
                           "time.sleep(60)"], start_new_session=True)
time.sleep(1.5)
group = os.getpgid(parent.pid)
L.register_child(parent)
L.stop_children()
time.sleep(1.0)
check("a process the app started is stopped with it", parent.poll() is not None, True)
try:
    os.killpg(group, 0)
    alive = True
except (ProcessLookupError, PermissionError):
    alive = False
check("and so is the process it started", alive, False)

print()
if fails:
    print(f"{len(fails)} failed:")
    for f in fails:
        print("  " + f)
    sys.exit(1)
print("all lifecycle checks passed")
