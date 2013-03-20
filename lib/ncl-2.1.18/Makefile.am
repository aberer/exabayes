SUBDIRS = ncl $(MAYBE_TEST_DIR) $(MAYBE_SWIG_DIR) #  example

pkgconfigdir= $(libdir)/pkgconfig
pkgconfig_DATA= nclv2.1.pc

EXTRA_DIST =  Jamroot  INSTALL COPYING AUTHORS NEWS BSDLicense.txt # data html test swig

ACLOCAL_AMFLAGS = -I config


dist-hook:
	find "$(distdir)/data" -depth -name .svn -and -type d -and -exec rm -rf {} \;
	find "$(distdir)/html" -depth -name .svn -and -type d -and -exec rm -rf {} \;
	find "$(distdir)/test" -depth -name .svn -and -type d -and -exec rm -rf {} \;
	find "$(distdir)/swig" -depth -name .svn -and -type d -and -exec rm -rf {} \;

