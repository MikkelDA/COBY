#!/bin/bash
# test_compat.sh - Run from inside github_version/
# Tests COBY import across multiple Python versions using micromamba

VERSIONS=(
    "3.9"
    "3.10"
    "3.11"
    "3.12"
    "3.13"
    "3.14"
)
PASS=()
FAIL=()

# Pass SRC_PATH as environment variable to avoid issues with spaces in path
export COBY_SRC="$(pwd)/src"

export MAMBA_ROOT_PREFIX="$HOME/micromamba"
eval "$(/home/au613219/.local/bin/micromamba shell hook --shell bash)"

for VERSION in "${VERSIONS[@]}"; do
    ENV_NAME="coby_test_py${VERSION//./}"
    echo ""
    echo "=============================="
    echo " Testing Python $VERSION"
    echo "=============================="

    # Create a minimal env if it doesn't exist
    if ! micromamba env list | grep -q "$ENV_NAME"; then
        echo "Creating environment $ENV_NAME..."

        micromamba create -y -n "$ENV_NAME" python="$VERSION" -c conda-forge
        if [ $? -ne 0 ]; then
            echo "SKIP: Python $VERSION could not be created"
            FAIL+=("Python $VERSION (env creation failed)")
            continue
        fi

        echo "Installing pip dependencies..."
        micromamba run -n "$ENV_NAME" pip install alphashape>=1.3.1 shapely>=2.0.0 matplotlib numpy scipy
        if [ $? -ne 0 ]; then
            echo "WARN: Some pip dependencies failed to install for Python $VERSION (may affect import test)"
        fi
    # else
        # echo "Environment $ENV_NAME already exists, skipping creation."
    fi

    # Try importing COBY from local src/ without installing it
    OUTPUT=$(micromamba run -n "$ENV_NAME" python -c "
import sys, os, compileall
sys.path.insert(0, os.environ['COBY_SRC'])

# First check all files compile cleanly
success = compileall.compile_dir(os.environ['COBY_SRC'], quiet=2, force=True)
if not success:
    print('SYNTAX ERROR: compilation failed')
    sys.exit(1)

try:
    import COBY
    print('PASS')
except SyntaxError as e:
    print(f'SYNTAX ERROR: {e.filename}:{e.lineno}: {e.msg}')
    print(f'  {e.text}')
except ImportError as e:
    print(f'IMPORT ERROR (likely missing dep): {e}')
except Exception as e:
    print(f'OTHER ERROR: {type(e).__name__}: {e}')
" 2>&1)

    echo "$OUTPUT"

    if echo "$OUTPUT" | grep -q "^PASS"; then
        PASS+=("Python $VERSION")
    else
        FAIL+=("Python $VERSION")
    fi
done

echo ""
echo "=============================="
echo " SUMMARY"
echo "=============================="
for v in "${PASS[@]}"; do echo "  ✓ $v"; done
for v in "${FAIL[@]}"; do echo "  ✗ $v"; done

# Exit with error code if any failures
if [ ${#FAIL[@]} -ne 0 ]; then
    exit 1
else
    exit 0
fi