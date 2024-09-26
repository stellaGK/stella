# See README.md for more information on these automated stella tests

# Test directory
test_dir = $(GK_HEAD_DIR)/tests

# Create a virtualenv and install the test dependencies there
create-test-virtualenv:
	bash -c "cd $(test_dir); \
		python3 -m venv venv; \
		source venv/bin/activate; \
		pip install -r requirements.txt; \
		deactivate; \
	"
	
# Check that all the Python dependencies are importable
check-automated-stella-tests-dependencies:
	@python3 $(test_dir)/check_automated_stella_tests_requirements.py
	
# Run the automated stella tests (the -r flag shows the full stella output and -rf shows extra info about the failed outputs)
run-automated-numerical-tests-for-stella: check-automated-stella-tests-dependencies $(GK_HEAD_DIR)/stella
	@echo "Running: automated stella tests"
	@pytest -v -rF -W ignore::DeprecationWarning
	
# Run the automated stella tests (the -r flag shows the full stella output)
# -rA shows all outputs, -rf shows only failed outputs 
run-automated-numerical-tests-for-stella-verbose: check-automated-stella-tests-dependencies $(GK_HEAD_DIR)/stella
	@echo "Running: automated stella tests printing all the output"
	@pytest -v -rA -W ignore::DeprecationWarning
