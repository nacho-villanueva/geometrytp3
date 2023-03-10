dnl Process this file with autoconf to produce a configure script.
dnl ------------------------------------------------------------------------
dnl initialisation
dnl ------------------------------------------------------------------------

dnl ./configure: sh internal 2K buffer overflow on HP-UX 9.xx
dnl thus, updating cache ./config.cache avoided.
define([AC_CACHE_LOAD], )dnl
define([AC_CACHE_SAVE], )dnl

AC_INIT(gmm, 5.4.2)
AC_CONFIG_MACRO_DIR([m4])
AC_CONFIG_HEADERS(config.h)
AX_PREFIX_CONFIG_H(src/gmm/gmm_arch_config.h,GMM)
AC_PREREQ(2.61)
AC_ARG_PROGRAM

PACKAGE="gmm"
MAJOR_VERSION="5"
MINOR_VERSION="4"
PATCH_VERSION="2"
AC_DEFINE_UNQUOTED([MAJOR_VERSION],5,[Major version number])
AC_DEFINE_UNQUOTED([MINOR_VERSION],4,[Minor version number])
AC_DEFINE_UNQUOTED([PATCH_VERSION],2,[Patch version number])
dnl VERSION=$MAJOR_VERSION.$MINOR_VERSION
VERSION=$MAJOR_VERSION.$MINOR_VERSION
echo "configuring $PACKAGE $VERSION..."

dnl ------------------------------------------------------------------------
dnl   init automake
dnl ------------------------------------------------------------------------

AM_INIT_AUTOMAKE([1.11 parallel-tests])

dnl -----------------------------------------------
dnl test du c++
dnl -----------------------------------------------

USER_CXXFLAGS="$CXXFLAGS"
AC_PROG_CXX(cxx KCC CC cc++ xlC aCC g++ c++ icc)
AC_PROG_CXXCPP
CXXFLAGS="${USER_CXXFLAGS}"
SUPLDFLAGS=""

AC_LANG_CPLUSPLUS

if test "x$prefix" = "xNONE"; then
  GFPREFIX=/usr/local;
else
  GFPREFIX="$prefix";
fi;

dnl AC_CXX_FULL_SPECIALIZATION_SYNTAX (c)Luc Maisonobe v 1.1.1.1 (2001/07/26)
dnl with some modification to test partial specialization
AC_CACHE_CHECK(whether the compiler recognizes the partial specialization syntax,
ac_cv_cxx_partial_specialization_syntax,
[AC_LANG_SAVE
 AC_LANG_CPLUSPLUS
 AC_TRY_COMPILE([
template<class T> class A        { public : int f () const { return 1; } };
template<class T> class A<T*>    { public:  int f () const { return 0; } };],[
A<float*> a; return a.f();],
 ac_cv_cxx_partial_specialization_syntax=yes, ac_cv_cxx_partial_specialization_s
yntax=no)
 AC_LANG_RESTORE
])
if test "$ac_cv_cxx_partial_specialization_syntax" != yes; then
  echo "Your compiler ($CXX) does not support partial template specialization, trash it"
  exit 1;
fi

AC_CANONICAL_HOST

echo "you are compiling gmm on a $host"

case $CXX in
 cxx)
        echo "Using Compaq cxx compiler"
        echo "WARNING : Control that you have at least Compaq C++ V6.3"
        here=`pwd`
        cd $srcdir
dnl     il faut utiliser -tweak au lieu des repositories ...
        CXXFLAGS="$CXXFLAGS -tweak -std strict_ansi -fast -Wl,-S -nopure_cname"
dnl     CXXFLAGS="$CXXFLAGS -ptr `pwd`/cxx_repository -std strict_ansi -O3"
        cd $here
        ;;
 CC)
        case $host in
        *irix*)
                echo "Using MIPSPRO CC on IRIX  (LD is set to CC)"
                LD=CC
dnl             CXXFLAGS="$CXXFLAGS -LANG:std -O3 -OPT:Olimit=0:roundoff=3:div_split=ON:alias=typed -TARG:platform=ip25"
                CXXFLAGS="$CXXFLAGS -LANG:std  -O3 "
dnl             CXXFLAGS="$CXXFLAGS -LANG:std  -O3 -ansiW "
                SUPLDFLAGS="-lCio"
                ;;
        *sun*)
                echo "Using SUN C++ WorkShop Compiler"
                CXXFLAGS="$CXXFLAGS +w2 -O3"
                ;;
        esac
        ;;
 aCC)
        echo "Using HP ANSI C++ Compiler aCC"
        CXXFLAGS="$CXXFLAGS -AA -fast"
        ;;
 g++* | c++)
        GCCVER=`$CXX --version | head -1 | cut -d ' ' -f3`
        echo "Using the GNU g++ compiler $GCCVER"
        case $GCCVER in
          2.95*)
                WSHADOW=""
                ;;
          *)
                WSHADOW="-Wshadow"
                ;;
        esac
        CXXFLAGS="$CXXFLAGS -ftemplate-depth-40 -pedantic -O3 -Wall -W $WSHADOW -Wpointer-arith -Wcast-qual -Wwrite-strings -Wconversion -Wredundant-decls -Wno-long-long"
        ;;
 icc | icpc)
        echo "Using INTEL icc"
dnl -tpp6 is for pentiumII and more
dnl -Xc is for ansi conformance
        CXXFLAGS="$CXXFLAGS -O3 -tpp6"
        ;;
 *)
        echo "Using a unknown compiler"
        CXXFLAGS="$CXXFLAGS -O3"
        ;;
esac

AC_SUBST(SUPLDFLAGS)

dnl ---------------------------OPENMP------------------------------
useopenmp=0
AC_ARG_ENABLE(openmp,
   [AS_HELP_STRING([--enable-openmp],[enable the multithreaded version of GMM++])],
   [ case $enableval in
        yes | "") useopenmp=YES ;;
        no) useopenmp=NO ;;
        *) AC_MSG_ERROR([bad value ${enableval} for --enable-openmp]) ;;
     esac],
   [useopenmp=NO]
)
if test x$useopenmp = xYES; then
  AC_OPENMP
  if test "x$ac_cv_prog_cxx_openmp" != "xunsupported" && test "x$ac_cv_prog_cxx_openmp" != "x"; then
    AC_SUBST(AM_CXXFLAGS,"$OPENMP_CXXFLAGS")
    CPPFLAGS="$CPPFLAGS -DGMM_HAVE_OPENMP"
  else
    AC_MSG_ERROR([OpenMP support not found. Use --enable-openmp=no flag to compile GMM++ without OpenMP]);
  fi
fi;
dnl ---------------------------END OF OPENMP-----------------------
dnl ------------------------------------------------------------------------
dnl   init libtools for shared libraries
dnl ------------------------------------------------------------------------

dnl AC_DISABLE_FAST_INSTALL

AM_ENABLE_STATIC

dnl AM_PROG_LIBTOOL
AM_PROG_LIBTOOL
AC_SUBST([LIBTOOL_DEPS])

AC_CHECK_HEADERS(sys/times.h)

BUILDER=`whoami`
AC_SUBST(BUILDER)
BUILDDATE=`date +%D,%H:%M:%S`
AC_SUBST(BUILDDATE)
CONFIGURE_ARGS=$ac_configure_args
AC_SUBST(CONFIGURE_ARGS)
LIBTOOL_VERSION_INFO="-version-info ${MAJOR_VERSION}:${MINOR_VERSION}:0"
AC_SUBST(LIBTOOL_VERSION_INFO)

dnl AC_CHECK_PROGS(RANLIB, ranlib)

dnl -----------------------------------------------
dnl sorties
dnl -----------------------------------------------
AC_CONFIG_FILES(\
        Makefile \
        tests/Makefile \
        include/Makefile gmm-config)
AC_OUTPUT
chmod a+x gmm-config
chmod a+x gmm-config
