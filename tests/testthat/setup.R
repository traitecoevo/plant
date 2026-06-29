# Escalate partial argument matching from a warning to an error for the
# entire test suite.  This catches misspelled parameter names that R would
# otherwise silently resolve via partial matching.
options(warnPartialMatchArgs = TRUE, warn = 2)
