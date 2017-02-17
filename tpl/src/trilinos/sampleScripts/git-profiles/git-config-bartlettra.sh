# Set your username
git config --global user.name "Roscoe A. Bartlett"
# Set your email address
git config --global user.email rabartl@sandia.gov
# Use colorized output when it makes sense
git config --global color.ui true
# Push only the current tracking branch (upstream in newer git versions)
git config --global push.default tracking
# Enable rerere by default
git config --global rerere.enabled 1
# Set up some shortcut commands
_SCRIPT_DIR=`echo $0 | sed "s/\(.*\)\/.*\.sh/\1/g"`
$_SCRIPT_DIR/../../cmake/tribits/python_utils/git-config-alias.sh
