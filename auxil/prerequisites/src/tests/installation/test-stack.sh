# shellcheck shell=bash
functions_exist() {
  : "${1?'Missing function name list'}"

  for function_name; do
    if [[ "$(type -t "${function_name}")" != "function" ]]; then
      emergency "function ${function_name} does not exist"
    fi
  done
}

detect_missing_variable_name() {
  : "${1?'detect_missing_variable name: Missing function name'}"

  stack_new __dummy_stack

  for function_name; do
    expected_error="Missing variable name in ${function_name}"
    error_message="$(${function_name} _dummy_stack_name 2>&1 >/dev/null)" || true
    error_message="${error_message##*: }" # grab text after final colon/space sequence
    if [[ "${error_message}" != "${expected_error}" ]]; then
      emergency "${function_name} failed to detect '${expected_error}'"
    fi
  done
  stack_destroy __dummy_stack
}

detect_missing_stack_name() {
  : "${1?'detect_missing_name: Missing function name'}"

  expected_error="Missing stack name"
  debug "detect_missing_stack_name: \${expected_error}=${expected_error}"

  for function_name; do
    debug "detect_missing_stack_name: executing \"\$(${function_name} 2>&1 > /dev/null)\" || true"
    error_message="$(${function_name} 2>&1 >/dev/null)" || true
    error_message="${error_message##*: }" # grab text after final colon/space sequence
    debug "detect_missing_stack_name: \${error_message}=${error_message}"
    if [[ "${error_message}" != "${expected_error}" ]]; then
      emergency "${function_name} failed to detect '${expected_error}'"
    fi
  done
}

detect_no_such_stack() {
  : "${1?'detect_no_such_name: Missing function name'}"

  stack_name="__nonexistent_stack"
  expected_error="No such stack -- "
  debug "detect_no_such_stack: \${expected_error}=${expected_error}"

  # shellcheck disable=SC2034
  dummy_variable=""

  for function_name; do
    debug "detect_no_such_stack: executing \"\$(${function_name} ${stack_name} dummy_variable 2>&1 > /dev/null)\" || true"
    error_message="$(${function_name} ${stack_name} dummy_variable 2>&1 >/dev/null)" || true
    missing_stack="${error_message##${expected_error}}" # grab text after final colon/space sequence
    debug "detect_no_such_stack: \${missing_stack}=${missing_stack}"
    if [[ "${missing_stack}" != "${stack_name}" ]]; then
      emergency "${function_name} failed to detect missing '${stack_name}'"
    fi
  done
}

detect_duplicate_stack_creation() {
  : "${1?'detect_duplicate_stack_creation: Missing function name'}"

  stack_name="__dummy_stack"
  stack_new "${stack_name}"

  expected_error="Stack already exists -- "
  debug "detect_duplicate_stack_creation: \${expected_error}=${expected_error}"
  function_name="stack_new"

  debug "detect_duplicate_stack_creation: executing \"\$(${function_name} ${stack_name} 2>&1 > /dev/null)\" || true"
  error_message="$(${function_name} ${stack_name} dummy_variable 2>&1 >/dev/null)" || true
  duplicate_stack="${error_message##${expected_error}}" # grab text after final colon/space sequence
  debug "detect_no_such_stack: \${duplicate_stack}=${duplicate_stack}"
  if [[ "${duplicate_stack}" != "${stack_name}" ]]; then
    emergency "${function_name} failed to detect duplicate '${stack_name}'"
  fi
  stack_destroy ${stack_name}
}

verify_stack_size_changes() {
  stack_new foobar
  stack_size foobar __foobar_size
  expected_size=0
  # shellcheck disable=SC2154
  if [[ "${__foobar_size}" != "${expected_size}" ]]; then
    emergency "verify_stack_size_changes: size=${__foobar_size} (expected ${expected_size})"
  fi

  stack_push foobar kernel
  stack_size foobar __foobar_new_size
  ((expected_size += 1))
  # shellcheck disable=SC2154
  if [[ "${__foobar_new_size}" != "${expected_size}" ]]; then
    emergency "verify_stack_size_changes: size=${__foobar_new_size} (expected 1)"
  fi

  stack_peek foobar tmp
  # shellcheck disable=SC2154
  if [[ "${tmp}" != "kernel" ]]; then
    emergency "verify_stack_size_changes: peeked item ('${tmp}') mismatch with pushed item ('kernel')"
  fi

  stack_size foobar __should_be_unchanged
  # shellcheck disable=SC2154
  if [[ "${__should_be_unchanged}" != "${expected_size}" ]]; then
    emergency "verify_stack_size_changes: size=${__should_be_unchanged} (expected ${expected_size})"
  fi

  stack_pop foobar popped
  # shellcheck disable=SC2154
  if [[ "${popped}" != "kernel" ]]; then
    emergency "verify_stack_size_changes: popped item ('${popped}') mismatch with pushed item ('kernel')"
  fi
  ((expected_size -= 1)) || true

  stack_size foobar __final_size
  # shellcheck disable=SC2154
  if [[ "${__final_size}" != "${expected_size}" ]]; then
    emergency "verify_stack_size_changes: size=${__final_size} (expected ${expected_size})"
  fi
}

test_stack() {
  # Verify availability of expected functions:
  functions_exist \
    stack_new stack_exists stack_print stack_push stack_peek stack_pop stack_size stack_destroy

  # Verify that each named function detects missing stack-name arguments:
  detect_missing_stack_name \
    stack_new stack_exists stack_print stack_push stack_peek stack_pop stack_size stack_destroy

  # Verify that each named function detects missing names of the variable to push, pop, or peek:
  detect_missing_variable_name \
    stack_pop stack_push stack_peek

  # Verify that each named function detects non-existent stacks:
  detect_no_such_stack \
    stack_destroy stack_peek stack_print stack_pop stack_push stack_size

  # Verify that duplicate creation generates the expected error:
  detect_duplicate_stack_creation \
    stack_new

  # Verify that push, peek, and pop yield correct size changes or lack thereof:
  verify_stack_size_changes

  info "test-stack.sh: All tests passed."

}
