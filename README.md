# icFlow3
tested with Qt 5.15.2 on Ubuntu 21.04

## Required libraries:

libgomp1

libtbb-dev

libhdf5-dev

libgmsh-dev

libeigen3-dev

libboost-dev

libboost-serialization-dev

libvtk9-dev

libvtk9-qt-dev

Additionally, Intel MKL libraries must be installed and the path should be updated manually in CMakeLists.txt.

If building from Qt Creator, ensure that the following Build Settings are configured consistently (so that Qt libraries have the same version and location): Qt5Charts_DIR, Qt5Core_DIR, Qt5Gui_DIR, Qt5Widgets_DIR, and Qt5_DIR.

![icyFlow3 screenshot](/screenshot.png?raw=true)

