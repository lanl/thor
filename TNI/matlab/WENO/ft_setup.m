%
fdwenolib.basedir = pwd;
%
fdwenolib.srcdir = fdwenolib.basedir + "/src";
%
fdwenolib.tt    = fdwenolib.srcdir + "/tt";
fdwenolib.ft    = fdwenolib.srcdir + "/ft";
fdwenolib.misc  = fdwenolib.srcdir + "/misc";
fdwenolib.ttlib = fileparts(fdwenolib.basedir)+"/external-libs/TT-toolbox";
%
addpath(fdwenolib.ft);
addpath(fdwenolib.misc);