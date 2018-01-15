Name:         he3lib
Version:      1.0
Release:      alt1

Summary:      He3 calculator, C/F/matlab/octave/cmdline interfaces
Group:        Sciences
License:      GPL
Url:          http://slazav.github.io/he3lib/index.html
Packager:     Vladislav Zavjalov <slazav@altlinux.org>

Source:       %name-%version.tar

BuildRequires: octave-devel libGraphicsMagick-c++-devel libGraphicsMagick-devel libgl2ps-devel libsuitesparse-devel
Requires: octave

%description
He3 calculator, C/F/matlab/octave/cmdline interfaces


%package devel
Group: Development/C
Summary: Development files of %name
Requires: %name = %version-%release

%description devel
This package contains the C headers and documentation required for building
programs based on %name.


%prep
%setup -q

%build
%makeinstall
%makeinstall install_octave

%files
%_bindir/he3
%_libdir/*.so
%dir %_datadir/octave/packages/he3lib
%dir %_datadir/octave/packages/he3lib/packinfo
%_datadir/octave/packages/he3lib/*.oct
%_datadir/octave/packages/he3lib/packinfo/*

%files devel
%_libdir/*.a
%_includedir/*


%changelog

