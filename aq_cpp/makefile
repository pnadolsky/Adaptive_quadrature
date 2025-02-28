# Detect OS
ifeq ($(OS),Windows_NT)
    $(info Windows)
    ifeq ($(shell uname -o 2>NUL),Cygwin)
        $(info Cygwin)
        IS_CYGWIN = 1
    else
        $(info Not Cygwin)
        IS_CYGWIN = 0
    endif
else
    $(info Something Else)
    IS_CYGWIN = 0
endif

# Set commands based on OS
ifeq ($(IS_CYGWIN),1)  # Cygwin
    $(info set commands Cygwin)
    MKDIR_BUILD = mkdir -p $(BUILD_DIR)
    MKDIR_BIN = mkdir -p $(BIN_DIR)
    RM = rm -rf $(BUILD_DIR) $(BIN_DIR)
else ifeq ($(OS),Windows_NT)  # Windows CMD
    $(info Set Commands Windows_NT)
    MKDIR_BUILD = @if not exist $(BUILD_DIR) ( echo "Creating $(BUILD_DIR)" & mkdir $(BUILD_DIR) )
    MKDIR_BIN = @if not exist $(BIN_DIR) ( echo "Creating $(BIN_DIR)" & mkdir $(BIN_DIR) )
    RM = if exist $(BUILD_DIR) rmdir /S /Q $(BUILD_DIR) & if exist $(BIN_DIR) rmdir /S /Q $(BIN_DIR)
else  # Linux/macOS
    $(info set commands Linux/macOS)
    MKDIR_BUILD = mkdir -p $(BUILD_DIR)
    MKDIR_BIN = mkdir -p $(BIN_DIR)
    RM = rm -rf $(BUILD_DIR) $(BIN_DIR)
endif

# Compiler and flags
CXX = g++
CXXFLAGS = -Iinclude -std=c++17 -Wall -Wextra -O2

# Directories
SRC_DIR = source
TEST_DIR = test
BUILD_DIR = build
BIN_DIR = bin

# Source and test files
SRC_FILES = $(wildcard $(SRC_DIR)/*.cpp)
TEST_FILES = $(wildcard $(TEST_DIR)/*.cpp)

# Object files
OBJ_FILES = $(patsubst $(SRC_DIR)/%.cpp, $(BUILD_DIR)/%.o, $(SRC_FILES))
TEST_OBJ_FILES = $(patsubst $(TEST_DIR)/%.cpp, $(BUILD_DIR)/%.o, $(TEST_FILES))

# Executable test targets
TEST_EXECUTABLES = $(patsubst $(TEST_DIR)/%.cpp, $(BIN_DIR)/%, $(TEST_FILES))

# Default target: build all tests
all: $(TEST_EXECUTABLES)

# Compile source files
$(BUILD_DIR)/%.o: $(SRC_DIR)/%.cpp
	$(MKDIR_BUILD)
	$(CXX) $(CXXFLAGS) -c $< -o $@

# Compile test files
$(BUILD_DIR)/%.o: $(TEST_DIR)/%.cpp
	$(MKDIR_BUILD)
	$(CXX) $(CXXFLAGS) -c $< -o $@

# Link test executables
$(BIN_DIR)/%: $(BUILD_DIR)/%.o $(OBJ_FILES)
	$(MKDIR_BIN)
	$(CXX) $(CXXFLAGS) $^ -o $@

# Clean build artifacts
clean:
	$(RM)

# Run all tests
test: all
	@for test in $(TEST_EXECUTABLES); do \
		echo "Running $$test"; \
		./$$test; \
	done

.PHONY: all clean test
