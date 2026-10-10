# config.sh - read by run_chain_v2.sh (BEAST X, settings A and B, round of October 2026)
# directory of the launcher 'beast' of the conda environment 'beast' (the launcher starts Java with -Xms64m -Xmx2048m)
BEAST_BIN=PATH_OF_THE_BEAST_X_ENVIRONMENT/bin
# as in the run script of the earlier round (work/beast/run_beast_chains.sh): limits the threads of the
# garbage collector of Java; it does not change the chain
JAVA_TOOL_OPTIONS_VALUE="-XX:ParallelGCThreads=2 -XX:ConcGCThreads=1"
SAVE_EVERY=500000
POLL_SECONDS=30
