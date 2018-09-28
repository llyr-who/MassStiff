# MassStiff FEM Solver

For my masters program we had to complete a handfull of so-called "special topics". We were allocated roughly two weeks
worth of time to complete these special topics (and write up!). For one of my special topics I decided to build a FEM package.

As time was limited there are certain choices that were made that were no great. For example, I knew nothing about
computational geometry going into this special topic and since O'Rouke places a lot of time into Ear-Clipping I jumped
in with Ear-Clipping. (TERRIBLE DECISION!)

Regardless, it was a lot of work to develop something like this in such a short time frame. If the interest is there I
would like to work on this further (writing new FEM packages), so please feel free to contact me.

As I was short on space there is an accociated YouTube channel that was linked to help show the results

https://www.youtube.com/channel/UCOSpHsqVj2jfoL8hdNlRYrg

anthonygoddard5678 AT hotmail -- DOT -- com

In order to get this working in Fedora, we need to run two commands.

    First we need to run:
        yum install /usr/include/X11/xlib.h
        yum install libXtst-devel
        yum install libXrandr-devel
        yum install mesa-libGL-devel
        yum groupinstall "X Software Development"
