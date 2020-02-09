#!/usr/bin/env python
try:
    from setuptools import setup
    from setuptools.command.build_ext import build_ext
    build_py = None
except ImportError:
    from distutils.core import setup
    try:
        from distutils.command.build_py import build_py_2to3 as build_py
    except ImportError:
        from distutils.command.build_py import build_py

import atexit
import distutils.ccompiler
import glob
import os
import shutil
import subprocess
import sys
import tarfile
import tempfile
import versioneer

try:
    from urllib import urlretrieve
    from urllib2 import Request, urlopen, URLError
except:
    # Maybe Python 3?
    from urllib.request import Request, urlopen, urlretrieve
    from urllib.error import URLError

from distutils.core import Extension
from distutils.util import get_platform

if os.name == 'nt':
  import msvcrt
  import time
else:
  from select import select

from subprocess import Popen, PIPE
from textwrap import dedent

# Global version number. Keep the format of the next line intact.
VERSION = '0.7.1'

# Check Python's version info and exit early if it is too old
if sys.version_info < (2, 5):
    print("This module requires Python >= 2.5")
    sys.exit(0)

###########################################################################

LIBIGRAPH_FALLBACK_INCLUDE_DIRS = ['/usr/include/igraph', '/usr/local/include/igraph']
LIBIGRAPH_FALLBACK_LIBRARIES = ['igraph']
LIBIGRAPH_FALLBACK_LIBRARY_DIRS = []

###########################################################################

def cleanup_tmpdir(dirname):
    """Removes the given temporary directory if it exists."""
    if dirname is not None and os.path.exists(dirname):
        shutil.rmtree(dirname)

def create_dir_unless_exists(*args):
    """Creates a directory unless it exists already."""
    path = os.path.join(*args)
    if not os.path.isdir(path):
        os.makedirs(path)

def ensure_dir_does_not_exist(*args):
    """Ensures that the given directory does not exist."""
    path = os.path.join(*args)
    if os.path.isdir(path):
        shutil.rmtree(path)

def find_static_library(library_name, library_path):
    """Given the raw name of a library in `library_name`, tries to find a
    static library with this name in the given `library_path`. `library_path`
    is automatically extended with common library directories on Linux and Mac
    OS X."""

    variants = ["lib{0}.a", "{0}.a", "{0}.lib", "lib{0}.lib"]
    if is_unix_like():
        extra_libdirs = ["/usr/local/lib64", "/usr/local/lib",
                "/usr/lib64", "/usr/lib", "/lib64", "/lib"]
    else:
        extra_libdirs = []

    for path in extra_libdirs:
        if path not in library_path and os.path.isdir(path):
            library_path.append(path)

    for path in library_path:
        for variant in variants:
            full_path = os.path.join(path, variant.format(library_name))
            if os.path.isfile(full_path):
                return full_path

def first(iterable):
    """Returns the first element from the given iterable."""
    for item in iterable:
        return item
    raise ValueError("iterable is empty")

def get_output(command):
    """Returns the output of a command returning a single line of output."""
    p = Popen(command, shell=True, stdin=PIPE, stdout=PIPE, stderr=PIPE)
    p.stdin.close()
    p.stderr.close()
    line=p.stdout.readline().strip()
    p.wait()
    if type(line).__name__ == "bytes":
        line = str(line, encoding="utf-8")
    return line, p.returncode

def http_url_exists(url):
    """Returns whether the given HTTP URL 'exists' in the sense that it is returning
    an HTTP error code or not. A URL is considered to exist if it does not return
    an HTTP error code."""
    class HEADRequest(Request):
        def get_method(self):
            return "HEAD"
    try:
        response = urlopen(HEADRequest(url))
        return True
    except URLError:
        return False

def is_unix_like(platform=None):
    """Returns whether the given platform is a Unix-like platform with the usual
    Unix filesystem. When the parameter is omitted, it defaults to ``sys.platform``
    """
    platform = platform or sys.platform
    platform = platform.lower()
    return platform.startswith("linux") or platform.startswith("darwin") or \
            platform.startswith("cygwin")

def preprocess_fallback_config():
    """Preprocesses the fallback include and library paths depending on the
    platform."""
    global LIBIGRAPH_FALLBACK_INCLUDE_DIRS
    global LIBIGRAPH_FALLBACK_LIBRARY_DIRS
    global LIBIGRAPH_FALLBACK_LIBRARIES
    
    if os.name == 'nt' and distutils.ccompiler.get_default_compiler() == 'msvc':
        # if this setup is run in the source checkout *and* the igraph msvc was build,
        # this code adds the right library and include dir
        version = '';
        if sys.version_info >= (2, 7) and sys.version_info < (3, 0):
          version = '27';
        elif sys.version_info >= (3, 4) and sys.version_info < (3, 5):
          version ='34';
        elif sys.version_info >= (3, 5):
          version ='35';
        all_msvc_dirs = glob.glob(os.path.join('..', 'igraph', 'igraph-*-msvc-py{0}'.format(version)))
        if len(all_msvc_dirs) > 0:
            if len(all_msvc_dirs) > 1:
                print("More than one MSVC build directory (igraph-*-msvc-py{0}) found!".format(version))
                print("It could happen that setup.py uses the wrong one! Please remove all but the right one!\n\n")

            msvc_builddir = all_msvc_dirs[-1]
            if not os.path.exists(os.path.join(msvc_builddir, "Release")):
                print("There is no 'Release' dir in the MSVC build directory\n(%s)" % msvc_builddir)
                print("Please build the MSVC build first!\n")
            else:
                print("Using MSVC build dir as a fallback: %s\n\n" % msvc_builddir)

                is_64bits = sys.maxsize > 2**32
                LIBIGRAPH_FALLBACK_INCLUDE_DIRS = [os.path.join(msvc_builddir, "include")]
                LIBIGRAPH_FALLBACK_LIBRARY_DIRS = [os.path.join(msvc_builddir, "Release")]

def version_variants(version):
    """Given an igraph version number, returns a list of possible version
    number variants to try when looking for a suitable nightly build of the
    C core to download from igraph.org."""

    result = [version]

    # Add trailing ".0" as needed to ensure that we have at least
    # major.minor.patch
    parts = version.split(".")
    while len(parts) < 3:
        parts.append("0")
        result.append(".".join(parts))

    return result

###########################################################################

class IgraphCCoreBuilder(object):
    """Class responsible for downloading and building the C core of igraph
    if it is not installed yet."""

    def __init__(self, versions_to_try, remote_url=None,
            show_progress_bar=True):
        self.versions_to_try = versions_to_try
        self.remote_url = remote_url
        self.show_progress_bar = show_progress_bar
        self._tmpdir = None

    @property
    def tmpdir(self):
        """The temporary directory in which igraph is downloaded and extracted."""
        if self._tmpdir is None:
            self._tmpdir = tempfile.mkdtemp(prefix="igraph.")
            atexit.register(cleanup_tmpdir, self._tmpdir)
        return self._tmpdir

    def download_and_compile(self):
        """Downloads and compiles the C core of igraph."""

        def _progress_hook(count, block_size, total_size):
            if total_size < 0:
                sys.stdout.write("\rDownloading %s... please wait." % local_file)
            else:
                percentage = count * block_size * 100.0 / total_size
                percentage = min(percentage, 100.0)
                sys.stdout.write("\rDownloading %s... %.2f%%" % (local_file, percentage))
            sys.stdout.flush()

        # Determine the remote URL if needed
        if self.remote_url is None:
            self.version, remote_url = self.find_first_version()
            if not self.version:
                print("Version %s of the C core of igraph is not found among the "
                        "nightly builds." % self.versions_to_try[0])
                print("Use the --c-core-version switch to try a different version.")
                print("")
                return False
            local_file = "igraph-%s.tar.gz" % self.version
        else:
            remote_url = self.remote_url
            local_file = remote_url.rsplit("/", 1)[1]

        # Now determine the full path where the C core will be downloaded
        local_file_full_path = os.path.join(self.tmpdir, local_file)

        # Download the C core
        if self.show_progress_bar:
            urlretrieve(remote_url, local_file_full_path, reporthook=_progress_hook)
            print("")
        else:
            print("Downloading %s... " % local_file)
            urlretrieve(remote_url, local_file_full_path)

        # Extract it in the temporary directory
        print("Extracting %s..." % local_file)
        archive = tarfile.open(local_file_full_path, "r:gz")
        archive.extractall(self.tmpdir)

        # Determine the name of the build directory
        self.builddir = None
        for name in os.listdir(self.tmpdir):
            full_path = os.path.join(self.tmpdir, name)
            if name.startswith("igraph") and os.path.isdir(full_path):
                self.builddir = full_path
                break

        if not self.builddir:
            print("Downloaded tarball did not contain a directory whose name "
                    "started with igraph; giving up build.")
            return False

        # Try to compile
        cwd = os.getcwd()
        try:
            print("Configuring igraph...")
            os.chdir(self.builddir)
            retcode = subprocess.call("CFLAGS=-fPIC CXXFLAGS=-fPIC ./configure --disable-tls --disable-gmp",
                    shell=True)
            if retcode:
                return False

            retcode = subprocess.call("make", shell=True)
            if retcode:
                return False

            libraries = []
            for line in open(os.path.join(self.builddir, "igraph.pc")):
                if line.startswith("Libs: ") or line.startswith("Libs.private: "):
                    words = line.strip().split()
                    libraries.extend(word[2:] for word in words if word.startswith("-l"))

            if not libraries:
                # Educated guess
                libraries = ["igraph"]

        finally:
            os.chdir(cwd)

        # Compilation succeeded; copy everything into igraphcore
        create_dir_unless_exists("igraphcore")
        ensure_dir_does_not_exist("igraphcore", "include")
        ensure_dir_does_not_exist("igraphcore", "lib")
        shutil.copytree(os.path.join(self.builddir, "include"),
                os.path.join("igraphcore", "include"))
        shutil.copytree(os.path.join(self.builddir, "src", ".libs"),
                os.path.join("igraphcore", "lib"))
        f = open(os.path.join("igraphcore", "build.cfg"), "w")
        f.write(repr(libraries))
        f.close()

        return True

    def find_first_version(self):
        """Finds the first version of igraph that exists in the nightly build
        repo from the version numbers provided in ``self.versions_to_try``."""
        for version in self.versions_to_try:
            remote_url = self.get_download_url(version=version)
            if http_url_exists(remote_url):
                return version, remote_url
        return None, None

    def get_download_url(self, version):
        return "http://igraph.org/nightly/get/c/igraph-%s.tar.gz" % version

    def run(self):
        return self.download_and_compile()


class BuildConfiguration(object):
    def __init__(self):
        global VERSION
        self.c_core_versions = version_variants(VERSION)
        self.c_core_url = None
        self.include_dirs = []
        self.library_dirs = []
        self.libraries = []
        self.extra_compile_args = []
        self.extra_link_args = []
        self.extra_objects = []
        self.show_progress_bar = True
        self.static_extension = False
        self.download_igraph_if_needed = True
        self.use_pkgconfig = True
        self._has_pkgconfig = None
        self.wait = True
        
    @property
    def has_pkgconfig(self):
        """Returns whether ``pkg-config`` is available on the current system
        and it knows about igraph or not."""
        if self._has_pkgconfig is None:
            if self.use_pkgconfig:
                line, exit_code = get_output("pkg-config igraph")
                self._has_pkgconfig = (exit_code == 0)
            else:
                self._has_pkgconfig = False
        return self._has_pkgconfig

    @property
    def build_ext(self):
        """Returns a class that can be used as a replacement for the
        ``build_ext`` command in ``distutils`` and that will download and
        compile the C core of igraph if needed."""
        try:
            from setuptools.command.build_ext import build_ext
        except ImportError:
            from distutils.command.build_ext import build_ext

        buildcfg = self
        class custom_build_ext(build_ext):
            def run(self):
                # Print a warning if pkg-config is not available or does not know about igraph
                if buildcfg.use_pkgconfig:
                    detected = buildcfg.detect_from_pkgconfig()
                else:
                    detected = False

                # Check whether we have already compiled igraph in a previous run.
                # If so, it should be found in igraphcore/include and
                # igraphcore/lib
                if os.path.exists("igraphcore"):
                    buildcfg.use_built_igraph()
                    detected = True

                # Download and compile igraph if the user did not disable it and
                # we do not know the libraries from pkg-config yet
                if not detected:
                    if buildcfg.download_igraph_if_needed and is_unix_like():
                        detected = buildcfg.download_and_compile_igraph()
                        if detected:
                            buildcfg.use_built_igraph()

                # Fall back to an educated guess if everything else failed
                if not detected:
                    buildcfg.use_educated_guess()

                # Replaces library names with full paths to static libraries
                # where possible
                if buildcfg.static_extension:
                    buildcfg.replace_static_libraries(exclusions=["m"])

                # Prints basic build information
                buildcfg.print_build_info()

                ext = first(extension for extension in self.extensions
                        if extension.name == "louvain._c_louvain")
                buildcfg.configure(ext)

                # Run the original build_ext command
                build_ext.run(self)

        return custom_build_ext

    def configure(self, ext):
        """Configures the given Extension object using this build configuration."""
        ext.include_dirs += self.include_dirs
        ext.library_dirs += self.library_dirs
        ext.libraries += self.libraries
        ext.extra_compile_args += self.extra_compile_args
        ext.extra_link_args += self.extra_link_args
        ext.extra_objects += self.extra_objects

    def detect_from_pkgconfig(self):
        """Detects the igraph include directory, library directory and the
        list of libraries to link to using ``pkg-config``."""
        if not buildcfg.has_pkgconfig:
            print("Cannot find the C core of igraph on this system using pkg-config.")
            return False

        cmd = "pkg-config igraph --cflags --libs"
        if self.static_extension:
            cmd += " --static"
        line, exit_code = get_output(cmd)
        if exit_code > 0 or len(line) == 0:
            return False

        opts = line.strip().split()
        self.libraries = [opt[2:] for opt in opts if opt.startswith("-l")]
        self.library_dirs = [opt[2:] for opt in opts if opt.startswith("-L")]
        self.include_dirs = [opt[2:] for opt in opts if opt.startswith("-I")]
        return True

    def download_and_compile_igraph(self):
        """Downloads and compiles the C core of igraph."""
        print("We will now try to download and compile the C core from scratch.")
        print("Version number of the C core: %s" % self.c_core_versions[0])
        if len(self.c_core_versions) > 1:
            print("We will also try: %s" % ", ".join(self.c_core_versions[1:]))
        print("")

        igraph_builder = IgraphCCoreBuilder(self.c_core_versions, self.c_core_url,
                show_progress_bar=self.show_progress_bar)
        if not igraph_builder.run():
            print("Could not download and compile the C core of igraph.")
            print("")
            return False
        else:
            return True

    def print_build_info(self):
        """Prints the include and library path being used for debugging purposes."""
        if self.static_extension:
            build_type = "static extension"
        else:
            build_type = "dynamic extension"
        print("Build type: %s" % build_type)
        print("Include path: %s" % " ".join(self.include_dirs))
        print("Library path: %s" % " ".join(self.library_dirs))
        print("Linked dynamic libraries: %s" % " ".join(self.libraries))
        print("Linked static libraries: %s" % " ".join(self.extra_objects))
        print("Extra compiler options: %s" % " ".join(self.extra_compile_args))
        print("Extra linker options: %s" % " ".join(self.extra_link_args))

    def process_args_from_command_line(self):
        """Preprocesses the command line options before they are passed to
        setup.py and sets up the build configuration."""
        # Yes, this is ugly, but we don't want to interfere with setup.py's own
        # option handling
        opts_to_remove = []
        for idx, option in enumerate(sys.argv):
            if not option.startswith("--"):
                continue
            if option == "--static":
                opts_to_remove.append(idx)
                self.static_extension = True
            elif option == "--no-download":
                opts_to_remove.append(idx)
                self.download_igraph_if_needed = False
            elif option == "--no-pkg-config":
                opts_to_remove.append(idx)
                self.use_pkgconfig = False
            elif option == "--no-progress-bar":
                opts_to_remove.append(idx)
                self.show_progress_bar = False
            elif option == "--no-wait":
                opts_to_remove.append(idx)
                self.wait = False                
            elif option.startswith("--c-core-version"):
                opts_to_remove.append(idx)
                if option == "--c-core-version":
                    value = sys.argv[idx+1]
                    opts_to_remove.append(idx+1)
                else:
                    value = option.split("=", 1)[1]
                self.c_core_versions = [value]
            elif option.startswith("--c-core-url"):
                opts_to_remove.append(idx)
                if option == "--c-core-url":
                    value = sys.argv[idx+1]
                    opts_to_remove.append(idx+1)
                else:
                    value = option.split("=", 1)[1]
                self.c_core_url = value

        for idx in reversed(opts_to_remove):
            sys.argv[idx:(idx+1)] = []

    def replace_static_libraries(self, exclusions=None):
        """Replaces references to libraries with full paths to their static
        versions if the static version is to be found on the library path."""
        if "stdc++" not in self.libraries:
            self.libraries.append("stdc++")

        if exclusions is None:
            exclusions = []

        for library_name in set(self.libraries) - set(exclusions):
            static_lib = find_static_library(library_name, self.library_dirs)
            if static_lib:
                self.libraries.remove(library_name)
                self.extra_objects.append(static_lib)

    def use_built_igraph(self):
        """Assumes that igraph is built already in ``igraphcore`` and sets up
        the include and library paths and the library names accordingly."""
        buildcfg.include_dirs = [os.path.join("igraphcore", "include")]
        buildcfg.library_dirs = [os.path.join("igraphcore", "lib")]
        buildcfg.static_extension = True

        buildcfg_file = os.path.join("igraphcore", "build.cfg")
        if os.path.exists(buildcfg_file):
            buildcfg.libraries = eval(open(buildcfg_file).read())

    def use_educated_guess(self):
        """Tries to guess the proper library names, include and library paths
        if everything else failed."""

        preprocess_fallback_config()

        global LIBIGRAPH_FALLBACK_LIBRARIES
        global LIBIGRAPH_FALLBACK_INCLUDE_DIRS
        global LIBIGRAPH_FALLBACK_LIBRARY_DIRS

        print("WARNING: we were not able to detect where igraph is installed on")
        print("your machine (if it is installed at all). We will use the fallback")
        print("library and include paths hardcoded in setup.py and hope that the")
        print("C core of igraph is installed there.")
        print("")
        print("If the compilation fails and you are sure that igraph is installed")
        print("on your machine, adjust the following two variables in setup.py")
        print("accordingly and try again:")
        print("")
        print("- LIBIGRAPH_FALLBACK_INCLUDE_DIRS",
            LIBIGRAPH_FALLBACK_INCLUDE_DIRS)
        print("- LIBIGRAPH_FALLBACK_LIBRARY_DIRS",
            LIBIGRAPH_FALLBACK_LIBRARY_DIRS)
        print("")

        seconds_remaining = 10 if self.wait else 0
        while seconds_remaining > 0:
            if seconds_remaining > 1:
                plural = "s"
            else:
                plural = ""

            sys.stdout.write("\rContinuing in %2d second%s; press Enter to continue "
                    "immediately. " % (seconds_remaining, plural))
            sys.stdout.flush()

            if os.name == 'nt':
              if msvcrt.kbhit():
                  if msvcrt.getch() == b'\r': # not '\n'
                      break
              time.sleep(1)
            else:
              rlist, _, _ = select([sys.stdin], [], [], 1)
              if rlist:
                  sys.stdin.readline()
                  break

            seconds_remaining -= 1
        sys.stdout.write("\r" + " "*65 + "\r")

        self.libraries = LIBIGRAPH_FALLBACK_LIBRARIES[:]
        if self.static_extension:
            self.libraries.extend(["xml2", "z", "m", "stdc++"])
        self.include_dirs = LIBIGRAPH_FALLBACK_INCLUDE_DIRS[:]
        self.library_dirs = LIBIGRAPH_FALLBACK_LIBRARY_DIRS[:]


###########################################################################

# Process command line options
buildcfg = BuildConfiguration();
buildcfg.process_args_from_command_line();

louvain_ext = Extension('louvain._c_louvain',
                    sources = glob.glob(os.path.join('src', '*.cpp')),
                    include_dirs=['include']);

cmdclass = versioneer.get_cmdclass()
cmdclass.update(build_ext=buildcfg.build_ext)

options =  dict(
  name = 'louvain',
  version=versioneer.get_version(),
  description = 'Louvain is a general algorithm for methods of community detection in large networks.',
  long_description=
    """
 Louvain is a general algorithm for methods of community detection in large networks.

 Please refer to the `documentation <http://louvain-igraph.readthedocs.io/en/latest>`_
 for more details.

 The source code of this package is hosted at `GitHub <https://github.com/vtraag/louvain-igraph>`_.
 Issues and bug reports are welcome at https://github.com/vtraag/louvain-igraph/issues.
    """,
  license = 'GPLv3+',
  url = 'https://github.com/vtraag/louvain-igraph',

  author = 'V.A. Traag',
  author_email = 'vincent@traag.net',
  test_suite = 'tests',

  provides = ['louvain'],
  package_dir = {'louvain': 'src'},
  packages = ['louvain'],
  ext_modules = [louvain_ext],
  install_requires = ['python-igraph >= {0}.0'.format(VERSION)],
  classifiers=[
      'Development Status :: 4 - Beta',
      'Environment :: Console',
      'Intended Audience :: End Users/Desktop',
      'Intended Audience :: Developers',
      'Intended Audience :: Science/Research',
      'License :: OSI Approved :: GNU General Public License v3 or later (GPLv3+)',
      'Operating System :: MacOS :: MacOS X',
      'Operating System :: Microsoft :: Windows',
      'Operating System :: POSIX',
      'Programming Language :: Python',
      'Programming Language :: C++',
      'Topic :: Scientific/Engineering :: Mathematics',
      'Topic :: Scientific/Engineering :: Information Analysis',
      'Topic :: Sociology'
    ],
  cmdclass=cmdclass,
);

if "macosx" in get_platform() and "bdist_mpkg" in sys.argv:
    # OS X specific stuff to build the .mpkg installer
    options["data_files"] = [ \
            ('/usr/local/lib', [os.path.join('..', '..', 'fatbuild', 'libigraph.0.dylib')])
    ]

if sys.version_info > (3, 0):
    if build_py is None:
        options["use_2to3"] = True
    else:
        options["cmdclass"]["build_py"] = build_py

setup(**options);
