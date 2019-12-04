# MassStiff FEM Solver

For my master's program we had to complete a handful of so-called "special topics". We were allocated roughly two weeks
worth of time to complete these special topics (including write up, which is listed as cppspecialtopic.pdf). For one of my special topics I decided to build a FEM package.

As time was limited there are certain choices that were made that were not great. For example, I knew nothing about
computational geometry going into this special topic and since O'Rouke places a lot of time into Ear-Clipping I jumped
in with Ear-Clipping. (TERRIBLE DECISION!)

Regardless, it was a lot of work to develop something like this in such a short time frame. If the interest is there, I
would like to work on this further (writing new FEM packages), so please feel free to open an issue.

As I was short on space there is an accociated <a href="https://www.youtube.com/channel/UCOSpHsqVj2jfoL8hdNlRYrg
" target="_blank"> YouTube channel</a> that was linked to help show the results.

At this stage I am revisiting the project to refactor and introduce modern C++ features.

In order to get this working in Fedora:

        yum install /usr/include/X11/xlib.h
        yum install libXtst-devel
        yum install libXrandr-devel
        yum install mesa-libGL-devel
        yum groupinstall "X Software Development"

To get things going in Ubuntu: 

        sudo apt-get install libx11-dev
        sudo apt-get install libxrender-dev
        sudo apt-get install mesa-common-dev
        sudo apt-get install libglu1-mesa-dev
