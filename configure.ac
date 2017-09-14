AC_INIT(asprin, 1.0.0)

# Allow the user to specify the path for Python
AC_ARG_WITH(python, [ --with-python       specify the location of python ])

# Check for python -- after this, PATH_TO_PYTHON will either contain a valid
# path to a location of the python executable, or (in the case where python can
# be found in the $PATH) will be empty.
with_python=`echo "$with_python" | sed "s|~|${HOME}|g"` # expand tilde
if test "$with_python" == ""
then
  with_python="$PATH"
  AC_SUBST(PATH_TO_PYTHON, "")
else
  AC_SUBST(PATH_TO_PYTHON, $with_python)
fi
AC_CHECK_PROG(PYTHON_CHECK,python,yes,no,$with_python)
if test x"$PYTHON_CHECK" != x"yes" ; then
  AC_MSG_ERROR([Could not find python, configuration failed])
fi

# set language to cpp
AC_PROG_CXX
AC_LANG(C++)

# Check for python
AC_CHECK_PROG(PYTHON_CHECK,$with_python,yes)
if test x"$PYTHON_CHECK" != x"yes" ; then
  AC_MSG_ERROR([Could not find the python executable, configuration failed])
fi

# Generate files from templates
AC_OUTPUT(Makefile)
