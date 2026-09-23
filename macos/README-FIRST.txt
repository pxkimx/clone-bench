Clone Bench — how to start
==========================

1. Unzip, then move "Clone Bench.app" anywhere (Applications or Desktop).

2. First open only — macOS Gatekeeper will say the app is "damaged" or "cannot be verified"
   because it is not signed with a paid Apple developer certificate. Fix it once:
      open Terminal, type   xattr -cr    (with a space after it), drag the app into the
      Terminal window so its path is filled in, press Return.
   Then double-click the app. (Alternative: right-click the app -> Open -> Open.)

3. The first launch installs its packages (Biopython, numpy, reportlab, the web server) into
   ~/Library/Application Support/CloneBench   (about a minute, needs internet once).
   Later launches take a few seconds and need no internet at all. Your browser opens
   http://localhost:8768 and the server keeps running quietly in the background. To stop it:
   "Quit Clone Bench" in the sidebar, or double-click "Stop Clone Bench.command".
   It uses its own port (8768), so it runs side by side with RNAseq Bench (8765) and
   MassSpec Bench (8766).

4. Everything is computed on this Mac. Only the optional assistant ("Ask Claude") talks to
   Anthropic's API; it reuses the API key you saved in RNAseq Bench, or paste one in Settings.

If the app will not open at all: double-click "Start Clone Bench.command" instead —
it runs the same thing inside a Terminal window so you can see every message.
The log is at ~/Library/Application Support/CloneBench/clone-bench.log

Needs Python 3.10–3.13 (python.org, Homebrew or Anaconda all work).
