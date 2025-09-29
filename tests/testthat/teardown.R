# Teardown code that runs after tests
if (exists("test_temp_dir") && dir.exists(test_temp_dir)) {
  unlink(test_temp_dir, recursive = TRUE)
}

# Reset any modified options
options(mestools.test_mode = NULL)
