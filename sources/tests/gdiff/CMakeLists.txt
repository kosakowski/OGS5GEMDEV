# This is Google diff-match-patch r103
# http://google-diff-match-patch.googlecode.com/

ADD_LIBRARY (GDiff gdiff.h gdiff.cpp)

SET_TARGET_PROPERTIES(GDiff PROPERTIES FOLDER Testing)