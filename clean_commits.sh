#!/bin/bash

# Script to remove Claude mentions from commit history

echo "Creating backup branch..."
git branch backup-before-cleanup-$(date +%Y%m%d_%H%M%S)

echo "Cleaning commit messages..."
git filter-branch --force --msg-filter '
    # Remove the robot emoji line and Claude co-author line
    sed -e "/🤖 Generated with/d" -e "/Co-Authored-By: Claude/d" | 
    # Remove any trailing empty lines
    sed -e :a -e "/^\s*$/{\$d;N;ba" -e "}"
' --tag-name-filter cat -- --all

echo "Cleanup complete!"
echo ""
echo "To push the cleaned history to GitHub (this will rewrite history!):"
echo "  git push --force-with-lease origin master"
echo ""
echo "To restore from backup if something went wrong:"
echo "  git reset --hard backup-before-cleanup-*"