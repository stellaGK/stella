# See README.md for more information on these automated stella tests
# The fortran tests have their own make file in the 'test_fortran_routines' subfolder

# Test directory
TEST_DIR = $(PROJECT_DIR)/AUTOMATIC_TESTS
.PHONY: create-test-virtualenv check-tests-dependencies numerical-tests numerical-tests-verbose benchmarks benchmarks-slow

# Create a virtualenv and install the test dependencies there
create-test-virtualenv:
	bash -c "cd $(TEST_DIR); \
		python3 -m venv venv; \
		source venv/bin/activate; \
		pip install -r requirements.txt; \
		deactivate; \
	"
	
# Check that all the Python dependencies are importable
check-tests-dependencies:
	@echo "Running: check-tests-dependencies"
	@python3 $(TEST_DIR)/check_automated_stella_tests_requirements.py
	
# Run the automated stella tests (the -r flag shows the full stella output and -rf shows extra info about the failed outputs)
numerical-tests: check-tests-dependencies $(PROJECT_DIR)/stella
	@echo "Running: automated stella tests"
	@python3 -m pytest -v -rF -W ignore::DeprecationWarning $(TEST_DIR)/numerical_tests/
	
# Run tests one by one
numerical-tests-1: check-tests-dependencies $(PROJECT_DIR)/stella
	@echo "Running: automated stella tests 1"
	@python3 -m pytest -v -rA -W ignore::DeprecationWarning $(TEST_DIR)/numerical_tests/test_1_whether_stella_runs
numerical-tests-2: check-tests-dependencies $(PROJECT_DIR)/stella
	@echo "Running: automated stella tests 2"
	@python3 -m pytest -v -rA -W ignore::DeprecationWarning $(TEST_DIR)/numerical_tests/test_2_geometry
numerical-tests-3: check-tests-dependencies $(PROJECT_DIR)/stella
	@echo "Running: automated stella tests 3"
	@python3 -m pytest -v -rA -W ignore::DeprecationWarning $(TEST_DIR)/numerical_tests/test_3_gyrokinetic_equation
numerical-tests-4: check-tests-dependencies $(PROJECT_DIR)/stella
	@echo "Running: automated stella tests 4"
	@python3 -m pytest -v -rA -W ignore::DeprecationWarning $(TEST_DIR)/numerical_tests/test_4_fluxtube
numerical-tests-5: check-tests-dependencies $(PROJECT_DIR)/stella
	@echo "Running: automated stella tests 5"
	@python3 -m pytest -v -rA -W ignore::DeprecationWarning $(TEST_DIR)/numerical_tests/test_5_diagnostics
numerical-tests-6: check-tests-dependencies $(PROJECT_DIR)/stella
	@echo "Running: automated stella tests 6"
	@python3 -m pytest -v -rF -W ignore::DeprecationWarning $(TEST_DIR)/numerical_tests/test_6_full_flux_surface
numerical-tests-7: check-tests-dependencies $(PROJECT_DIR)/stella
	@echo "Running: automated stella tests 7"
	@python3 -m pytest -v -rA -W ignore::DeprecationWarning $(TEST_DIR)/numerical_tests/test_7_electromagnetic
numerical-tests-system:
	@python3 $(TEST_DIR)/system_information.py
	
	
# Run the automated stella tests on older stella versions
# TODO-HT extend the other tests to the old stella releases
# @python3 -m pytest -v -rF -W ignore::DeprecationWarning --stella_version 0.5 $(TEST_DIR)/numerical_tests/test_2_geometry
numerical-tests-v0.5: check-tests-dependencies $(PROJECT_DIR)/stella
	@echo "Running: automated stella tests"
	@python3 -m pytest -v -rF -W ignore::DeprecationWarning --stella_version 0.5 $(TEST_DIR)/numerical_tests/test_1_whether_stella_runs
numerical-tests-v0.6: check-tests-dependencies $(PROJECT_DIR)/stella
	@echo "Running: automated stella tests"
	@python3 -m pytest -v -rF -W ignore::DeprecationWarning --stella_version 0.6 $(TEST_DIR)/numerical_tests/test_1_whether_stella_runs
numerical-tests-v0.7: check-tests-dependencies $(PROJECT_DIR)/stella
	@echo "Running: automated stella tests"
	@python3 -m pytest -v -rF -W ignore::DeprecationWarning --stella_version 0.7 $(TEST_DIR)/numerical_tests/test_1_whether_stella_runs
	
	
# Run the automated stella tests (the -r flag shows the full stella output and -rf shows extra info about the failed outputs)
pyrokinetics-tests: check-tests-dependencies $(PROJECT_DIR)/stella
	@echo "Running: automated stella tests for pyrokinetics"
	@python3 -m pytest -v -rF -W ignore::DeprecationWarning $(TEST_DIR)/pyrokinetics_tests/
	
# Run the automated stella tests (the -r flag shows the full stella output)
# -rA shows all outputs, -rf shows only failed outputs 
numerical-tests-verbose: check-tests-dependencies $(PROJECT_DIR)/stella
	@echo "Running: automated stella tests printing all the output"
	@python3 -m pytest -v -rA -W ignore::DeprecationWarning $(TEST_DIR)/numerical_tests
	
	# Run the automated stella tests (the -r flag shows the full stella output and -rf shows extra info about the failed outputs)
benchmarks: check-tests-dependencies $(PROJECT_DIR)/stella
	@echo "Running: automated stella benchmarks"
	@python3 -m pytest -v -rF -W ignore::DeprecationWarning $(TEST_DIR)/quick_benchmarks
	
# Run the automated stella tests (the -r flag shows the full stella output and -rf shows extra info about the failed outputs)
benchmarks-slow: check-tests-dependencies $(PROJECT_DIR)/stella
	@echo "Running: automated stella benchmarks"
	@python3 -m pytest -v -rF -W ignore::DeprecationWarning $(TEST_DIR)/slow_benchmarks
