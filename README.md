# skinObjectness
MATLAB code for data driven skin detection in cluttered USAR environments

The code provides an implementation of detecting skin regions in cluttered environments. Use the following tips to parse through the code.
** Update Yinitialize.m or directly put directory path in the skinObjectness.m to begin.
**The code is provided for generating skin probability maps and ROI windows for images. However, further implementation of the algorithm to video frames is also provided but needs some tweaking by the user.
** Use objectnessTracking.m for finding skin objectness in video frames by updating/enabling the 'directory' variable to parse through a directory of image frames in skinObjectness.m.
**The code uses some of the inbuilt functions in MATLAB and works fine in MATLAB with versions released after 2013.
