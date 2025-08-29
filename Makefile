# Compiler settings
CXX = g++
CXXFLAGS = -Wall -Werror -pedantic -std=gnu++98 -ggdb3

# Directories
SRC_DIR = src
CORE_DIR = $(SRC_DIR)/core
BUILD_DIR = build
INPUT_DIR = inputs

# Source files
CORE_SRCS = $(CORE_DIR)/mesh.cpp $(CORE_DIR)/triangle.cpp
MAIN_SRC = $(SRC_DIR)/main.cpp
SRCS = $(CORE_SRCS) $(MAIN_SRC)

# Object files (in build directory)
CORE_OBJS = $(BUILD_DIR)/mesh.o $(BUILD_DIR)/triangle.o
MAIN_OBJ = $(BUILD_DIR)/main.o
OBJS = $(CORE_OBJS) $(MAIN_OBJ)

# Target executable
TARGET = $(BUILD_DIR)/mesh-generator

# Default target
all: $(TARGET)

# Main executable
$(TARGET): $(OBJS) | $(BUILD_DIR)
	$(CXX) -o $@ $^

# Object file rules
$(BUILD_DIR)/main.o: $(MAIN_SRC) $(CORE_DIR)/mesh.hpp $(CORE_DIR)/node.hpp $(CORE_DIR)/triangle.hpp | $(BUILD_DIR)
	$(CXX) $(CXXFLAGS) -c $< -o $@

$(BUILD_DIR)/mesh.o: $(CORE_DIR)/mesh.cpp $(CORE_DIR)/mesh.hpp $(CORE_DIR)/triangle.hpp $(CORE_DIR)/node.hpp | $(BUILD_DIR)
	$(CXX) $(CXXFLAGS) -c $< -o $@

$(BUILD_DIR)/triangle.o: $(CORE_DIR)/triangle.cpp $(CORE_DIR)/triangle.hpp $(CORE_DIR)/node.hpp | $(BUILD_DIR)
	$(CXX) $(CXXFLAGS) -c $< -o $@

# Create build directory
$(BUILD_DIR):
	mkdir -p $(BUILD_DIR)

# Test targets
test: $(TARGET)
	@echo "Running test with simple-mesh.txt:"
	@./$(TARGET) --step1 $(INPUT_DIR)/simple-mesh.txt

test-all: $(TARGET)
	@echo "Testing all steps with simple-mesh.txt:"
	@echo "Step 1:"
	@./$(TARGET) --step1 $(INPUT_DIR)/simple-mesh.txt
	@echo -e "\nStep 2:"
	@./$(TARGET) --step2 $(INPUT_DIR)/simple-mesh.txt
	@echo -e "\nStep 3:"
	@./$(TARGET) --step3 $(INPUT_DIR)/simple-mesh.txt
	@echo -e "\nStep 4:"
	@./$(TARGET) --step4 $(INPUT_DIR)/simple-mesh.txt

# Install target (copy to root for backward compatibility)
install: $(TARGET)
	cp $(TARGET) ./mesh-generator

# Clean targets
clean:
	rm -rf $(BUILD_DIR)
	rm -f mesh-generator

clean-all: clean
	rm -f mesh-step* *.o *~

# Format code
format:
	clang-format -i $(SRC_DIR)/*.cpp $(CORE_DIR)/*.cpp $(CORE_DIR)/*.hpp

# Help target
help:
	@echo "Available targets:"
	@echo "  all        - Build the mesh generator (default)"
	@echo "  test       - Run basic test with step1"
	@echo "  test-all   - Run tests for all steps"
	@echo "  install    - Copy executable to root directory"
	@echo "  clean      - Remove build directory"
	@echo "  clean-all  - Remove all generated files"
	@echo "  format     - Format source code with clang-format"
	@echo "  help       - Show this help message"

.PHONY: all test test-all install clean clean-all format help