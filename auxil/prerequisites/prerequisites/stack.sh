# A stack, using bash arrays.
# ---------------------------------------------------------------------------
# This software is released under a BSD license, adapted from
# <http://opensource.org/licenses/bsd-license.php>
#
# Copyright &copy; 1989-2012 Brian M. Clapper.
# All rights reserved.
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are met:
#
# * Redistributions of source code must retain the above copyright notice,
#  this list of conditions and the following disclaimer.
#
# * Redistributions in binary form must reproduce the above copyright notice,
#  this list of conditions and the following disclaimer in the documentation
#  and/or other materials provided with the distribution.
#
# * Neither the name "clapper.org" nor the names of its contributors may be
#  used to endorse or promote products derived from this software without
#  specific prior written permission.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
# AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
# ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
# LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
# CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
# SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
# INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
# CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
# ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
# POSSIBILITY OF SUCH DAMAGE.
#

# Source: http://brizzled.clapper.org/blog/2011/10/28/a-bash-stack/

# Destroy a stack
#
# Usage: stack_destroy name
function stack_destroy
{
    : "${1?'Missing stack name'}"

    if no_such_stack "$1"
    then
        echo "No such stack -- $1" >&2
        return 1
    fi

    eval "unset _stack_$1 _stack_$1_i"
    return 0
}

# Push one or more items onto a stack.
#
# Usage: stack_push stack item ...
function stack_push
{
    : "${1?'Missing stack name'}"
    : "${2?'Missing variable name in stack_push'}"

    if no_such_stack "$1"
    then
        echo "No such stack -- $1" >&2
        return 1
    fi

    stack=$1
    shift 1

    while (( $# > 0 ))
    do
        eval '_i=$'"_stack_${stack}_i"
        eval "_stack_${stack}[$_i]='$1'"
        eval "let _stack_${stack}_i+=1"
        shift 1
    done

    unset _i
    return 0
}


# Get the size of a stack
#
# Usage: stack_size name var
#
# Example:
#    stack_size mystack n
#    echo "Size is $n"
function stack_size
{
    : "${1?'Missing stack name'}"
    : "${2?'Missing name of variable for stack size result'}"
    if no_such_stack "$1"
    then
        echo "No such stack -- $1" >&2
        return 1
    fi
    # TODO: revise the eval below to eliminate the need for this pop/push
    # sequene, which is a workaround to prevent an error that occurs with
    # if the stack is new and has not been the target of a stack_push.
    stack_push $1 __push_junk
    stack_pop $1 __pop_trash
    eval "$2"='$'"{#_stack_$1[*]}"
}

function no_such_stack
{
    : "${1?'Missing stack name'}"
    stack_exists "$1"
    ret=$?
    declare -i x
    let x="1-$ret" || true
    return $x
}

#### Functions modified by Damian Rouson

# These functions were modified to work with the shell settings in bash3boilerplate
# (https://github.com/zbeekman/bash3boilerplate), primarily the "set -o unset" and
# "set -o pipefail" settings.  The required modifications are described below.
#

# Pop the top element from the stack.
#
# Usage: stack_pop name var
#
# Example:
#    stack_pop mystack top
#    echo "Got $top"
#
# Modification for use with bash3boilerplate:
#    replaced "let _i-=1" with  "(( _i-=1 )) || true"

function stack_pop
{
    : "${1?'Missing stack name'}"
    : "${2?'Missing variable name in stack_pop'}"

    if no_such_stack "$1"
    then
        echo "No such stack -- $1" >&2
        return 1
    fi

    eval 'let _i=$'"_stack_$1_i"

    if [[ "$_i" -eq 0 ]]
    then
        echo "Empty stack -- $1" >&2
        return 1
    fi

   (( _i-=1 )) || true
    eval "$2"='$'"{_stack_$1[$_i]}"
    eval "unset _stack_$1[$_i]"
    eval "_stack_$1_i=$_i"
    unset _i
    return 0
}

# Print a stack to stdout.
#
# Usage: stack_print name
#
# Modification for use with bash3boilerplate:
#    replaced "let _i=${_i}-1" with "(( _i=${_i}-1 )) || true" to support execution with "set -o nounset"

function stack_print
{
    : "${1?'Missing stack name'}"

    if no_such_stack "$1"
    then
        echo "No such stack -- $1" >&2
        return 1
    fi

    tmp=""
    eval 'let _i=$'"_stack_$1_i"

    while (( _i > 0 ))
    do
        eval 'e=$'"{_stack_$1[$((--_i))]}" # pre-decrement
	# shellcheck  disable=SC2154
        tmp="$tmp $e"
    done
    # shellcheck  disable=SC2086
    echo "(" $tmp ")"
}

# Create a new stack.
#
# Usage: stack_new name
#
# Example: stack_new x
#
# Modification for use with bash3boilerplate:
#    added "|| true" to allow execution with "set -o nounset" on OS X

function stack_new
{
    : "${1?'Missing stack name'}"
    if stack_exists "$1"
    then
        echo "Stack already exists -- $1" >&2
        return 1
    fi

    eval "_stack_$1=()"
    eval "_stack_$1_i=0"

    variableName="_stack_$1_i"
    variableVal="0"
    eval "${variableName}"="$(echo -ne \""${variableVal}"\")"

    return 0
}

# Verify stack existence.
#
# Usage: stack_exists name
#
# Example: stack_new x
#
# Modification for use with bash3boilerplate:
#    added curly braces in eval statement to allow execution with "set -o nounset"

function stack_exists
{
    : "${1?'Missing stack name'}"

    eval '_i=$'"{_stack_$1_i:-}"
    if [[ -z "${_i:-}" ]]
    then
        return 1
    else
        return 0
    fi
}

#### Functions added by Damian Rouson

# Get the top element from the stack and return the stack
# to its state before invocation of the function.
#
# Usage: stack_peek name var
#
# Example:
#    stack_peek mystack top
#    echo "Got $top"
function stack_peek
{
  : "${1?'Missing stack name'}"
  : "${2?'Missing variable name in stack_peek'}"

  if no_such_stack "$1"
  then
      echo "No such stack -- $1" >&2
      return 1
  fi

  stack_pop "$1" "$2"
  eval argument_name="\$$2"
  # shellcheck disable=SC2154
  stack_push "$1" "$argument_name"
}
