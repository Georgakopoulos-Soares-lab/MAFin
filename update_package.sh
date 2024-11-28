#!/bin/bash

# Initialize Conda for non-interactive shells
CONDA_BASE=$(conda info --base)  # Get the base conda installation directory
source "$CONDA_BASE/etc/profile.d/conda.sh"  # Source Conda's shell integration

# Function to increment the minor version number
increment_version() {
    version=$1
    major=$(echo $version | cut -d. -f1)
    minor=$(echo $version | cut -d. -f2)
    patch=$(echo $version | cut -d. -f3)

    if [ -z "$patch" ]; then
        # If no patch version exists, treat as major.minor format
        new_minor=$((minor + 1))
        echo "$major.$new_minor"
    else
        # If patch version exists, treat as major.minor.patch format
        new_minor=$((minor + 1))
        echo "$major.$new_minor.$patch"
    fi
}

# Check if the Conda environment 'py38_deployment_pypi' exists using conda env list
env_exists=$(conda env list | grep -w py38_deployment_pypi)

# If environment doesn't exist, create it
if [ -z "$env_exists" ]; then
    echo "Conda environment 'py38_deployment_pypi' not found. Creating a new environment..."
    
    # Create the environment with Python 3.8 and install required packages
    conda create -n py38_deployment_pypi python=3.8 -y

    # Check if the environment creation succeeded
    if [ $? -ne 0 ]; then
        echo "Failed to create the Conda environment 'py38_deployment_pypi'. Exiting..."
        exit 1
    fi
fi

# Activate the Conda environment
# Workaround for conda activate in non-interactive shells
eval "$(conda shell.bash hook)"  # Properly initialize Conda's shell commands
conda activate py38_deployment_pypi

if [ $? -ne 0 ]; then
    echo "Failed to activate the Conda environment 'py38_deployment_pypi'. Exiting..."
    exit 1
fi

# Install required dependencies in case they are missing
conda install twine setuptools wheel -y

# Locate setup.py file in the current directory
if [ ! -f "setup.py" ]; then
    echo "setup.py not found in the current directory."
    exit 1
fi

# Extract the current version from setup.py
current_version=$(grep -oP "(?<=version=')[^']+" setup.py)
if [ -z "$current_version" ]; then
    echo "Could not find the version in setup.py"
    exit 1
fi

echo "Current version: $current_version"

# Increment the minor version
new_version=$(increment_version $current_version)

# Replace the version in setup.py with the incremented version
sed -i "s/version='$current_version'/version='$new_version'/" setup.py
echo "Updated version to $new_version"

# Build the package
python setup.py sdist bdist_wheel

# Upload to PyPI using Twine
twine upload dist/*

# Print completion message
echo "Package version $new_version has been successfully uploaded."
