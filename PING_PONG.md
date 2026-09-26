# Ping-Pong Compilation Protocol

## Mac Status (Written by Mac Claude)
```
STATUS: IDLE
LAST_ACTION: None
TIMESTAMP: 
MESSAGE: Waiting for Linux to compile
```

## Linux Status (Written by Linux Claude)
```
STATUS: IDLE
LAST_ACTION: None
TIMESTAMP:
ERRORS: None
MESSAGE: Ready to start
```

## Instructions

### For Linux Claude:
1. Compile all .c files in src/
2. Report errors to "Linux Status" section
3. Set STATUS: ERRORS or SUCCESS
4. Commit and push PING_PONG.md
5. Wait for Mac to fix

### For Mac Claude:
1. Monitor PING_PONG.md for Linux STATUS: ERRORS
2. Read error messages
3. Fix code, commit, push
4. Set Mac STATUS: FIXED
5. Wait for Linux to compile again

### Protocol:
- Linux: COMPILING → ERRORS or SUCCESS
- Mac: FIXING → FIXED → WAITING
- Repeat until SUCCESS
