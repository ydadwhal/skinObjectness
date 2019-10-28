function [common] = Yinitialize()
%Author: Yogeshwar Singh Dadwhal
%Date : 05/04/2016


curr=pwd;
%path of current directory

common.path=[curr '\'];

%current drive with \
common.currDir=curr(1,1:3);

%Directory for results
common.out=[common.currDir 'Results\'];

%Input Test Data Directory
% common.data='C:\Data\TestImages\'
% common.data='C:\Data\Clutter 2 flikr creative common liscense sleeping man\';
% common.data='C:\Data\Mohali01052016\';
% common.data='C:\Data\Mohali01052016\New folder\';
% common.data='C:\Data\SFA human skin image database\ORI\';
% common.data='C:\Data\Clutter 2 flikr creative common liscense sleeping man\';
common.data='C:\Users\Yogeshwar\Desktop\saliency data\Cheng and mitra human data\'
common.GT='F:\Data\SFA human skin image database\GT\';
%Input Skin Data

common.results=[curr '\Results\'];

common.skindata='H:\cdata\Data\skin segmentation dataset\';

addpath(genpath(curr));

addpath([common.currDir 'Dated MAT\feature\']);


