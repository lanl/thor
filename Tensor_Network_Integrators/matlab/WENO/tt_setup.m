%
fdwenolib.basedir = pwd;
%
fdwenolib.srcdir = fdwenolib.basedir + "/src";
%
fdwenolib.tt    = fdwenolib.srcdir + "/tt";
fdwenolib.ft    = fdwenolib.srcdir + "/ft";
fdwenolib.misc  = fdwenolib.srcdir + "/misc";
fdwenolib.ttlib = fileparts(fdwenolib.basedir)+"../utils/TT-toolbox";
%
addpath(fdwenolib.tt);
addpath(fdwenolib.misc);
%
run(fdwenolib.ttlib+"/setup.m");
