#!/bin/sh

WIKI_DIR=../wiki
FIGURE_DIR="figure"
GIT_WIKI="--git-dir=$WIKI_DIR/.git --work-tree=$WIKI_DIR"
SHORT_SHA=$(git rev-parse --short HEAD)

if test ! -d "$WIKI_DIR"
then
    WIKI_URI=$(git ls-remote --get-url | grep github.com: | head -n1 | sed 's/.git$/.wiki.git/')
    if test -z $WIKI_URI
    then
	echo "Error: could not determine github repo"
	exit 1
    else
	echo "bar"
    fi
    git clone $WIKI_URI ../wiki
fi

# Check that we don't have any residual "unnamed chunk" figures that I
# never want committed.
if [[ -d $FIGURE_DIR &&
      -n "$(find $FIGURE_DIR -maxdepth 1 -name '*unnamed*' -print -quit)" ]]
then
    echo "Error: Unnamed figure chunks found"
    exit 1
fi

# I'm not sure if this is ideal.  Could be over cautious, given we
# don't generally care all that much about the internal state of the
# wiki.
if test -n "$(git $GIT_WIKI status -s 2> /dev/null)"
then
    echo "Error: wiki git is dirty"
    echo "(consider commiting or 'git reset --hard HEAD' in the wiki dir)"
    exit 1
fi

while read S; do
    S_BASE="${S%.*}"

    echo "Updating ${S_BASE}"

    # Delete all old figures for this case.  This is needed because we
    # need to work out where images are no longer needed or they'll
    # persist in the wiki.  Perhaps better would be to use rsync with
    # some pattern matching, but we need to hit the 'git rm' at some
    # point.
    git $GIT_WIKI rm -f --quiet --ignore-unmatch -- "$FIGURE_DIR/${S_BASE}_*"

    # That *might* have deleted the wiki figure directory:
    mkdir -p $WIKI_DIR/$FIGURE_DIR

    # Copy the new stuff over:
    cp ${S_BASE}.md $WIKI_DIR
    cp $FIGURE_DIR/${S_BASE}_* $WIKI_DIR/$FIGURE_DIR/

    # And add it to git:
    git $GIT_WIKI add ${S_BASE}.md "$FIGURE_DIR/${S_BASE}_*"
done < .wiki_scripts

if git $GIT_WIKI status --porcelain --untracked-files=no | grep --quiet '^[A-Z]'
then
    git $GIT_WIKI commit -q -m "Updated wiki [at $SHORT_SHA]"
else
    echo "Wiki already up to date"
fi
