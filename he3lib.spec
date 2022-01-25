Name:         he3lib
Version:      1.0
Release:      alt1

Summary:      He3 calculator, C/F/matlab/octave/cmdline interfaces
Group:        Sciences/Physics
License:      GPL
Url:          http://slazav.github.io/he3lib/index.html
Packager:     Vladislav Zavjalov <slazav@altlinux.org>

Source:       %name-%version.tar

BuildRequires(pre): rpm-build-octave
BuildRequires: octave-devel gcc-fortran

#BuildRequires: octave-devel libGraphicsMagick-c++-devel libGraphicsMagick-devel
#BuildRequires: libgl2ps-devel libsuitesparse-devel libhdf5-devel libqrupdate-devel
#BuildRequires: libreadline-devel libpcre-devel
Requires: octave

%description
He3 calculator, C/F/matlab/octave/cmdline interfaces


%package devel
Group: Sciences/Physics
Summary: Development files of %name
Requires: %name = %version-%release

%description devel
This package contains the C headers and documentation required for building
programs based on %name.


%package -n octave-%name
Group: Sciences/Physics
Summary: Octave library of %name
Requires: %name = %version-%release
Requires: octave

%description -n octave-%name
This package contains %name Octave library.


%prep
%setup -q

%build
%make_build cmdline library octave-pkg

%install
%define _makeinstall_target install_cmdline install_library install_headers
%makeinstall
octave-cli --eval "pkg prefix %buildroot%octarchprefix;\
                   pkg install -nodeps -verbose -local octave-he3lib.tgz"

%files
%_bindir/he3
%_libdir/*.so

%files devel
%_includedir/*

%files -n octave-%name
%octarchprefix/%name-%version

%changelog

