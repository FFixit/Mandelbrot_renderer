# Mandelbrot Renderer
This is a project I did to improve my C/C++ skill after learning the basics from a book and also because the math is interesting.

![Demo Image](./demo_image.png "Mandelbrot Render")

It can render Bitmap images of the mandelbrot set, like the one above, with arbitrary zoom. There is also a function for creating a series of images that can be put together to make a video, like i did here: https://www.youtube.com/watch?v=MGPfdfrX7Bs

## How does this work?
This application can render images of the [mandelbrot set](https://en.wikipedia.org/wiki/Mandelbrot_set). Every pixel is mapped to a point in the complex plane. But instead of calculating the standard iteration (z_next = z²+c), I implemented two major optimisations:
-   [Perturbation](https://fractalwiki.org/wiki/Perturbation_theory): Instead of iterating for every pixel, only do the iteration for one reference point, and calculate every other pixel relative to the reference point.
-   [Series approximation](https://en.wikipedia.org/wiki/Plotting_algorithms_for_the_Mandelbrot_set#Perturbation_theory_and_series_approximation): Basically the math is transformed in a way, where some terms can be omitted, while the image is still the same (This is error-prone and my implementation isn't perfect)
## How to build?
This is a Visual Studio Project. To be able to build this, you must have the boost and GNU MPFR libraries installed and available to Visual Studio. To do this you have to:
-   Download '[boost](https://www.boost.org/)' and add the directory to `Solution Explorer > mandelbrot > Rightclick > Properties > Configuration Properties > C/C++ > Additional Include Directories`
-   Download or compile binaries for 'mpir' and 'mpfr' ([Can be downloaded here](http://www.holoborodko.com/pavel/wp-content/plugins/download-monitor/download.php?id=5))
    - Add the directories containing the headers to `Solution Explorer > mandelbrot > Rightclick > Properties > Configuration Properties > C/C++ > Additional Include Directories`
    - Add the directories containing the files `mpir.lib` and `mpfr.lib` to `Solution Explorer > mandelbrot > Rightclick > Properties > Configuration Properties > Linker > Additional Library Directories`
    - Depending on how you built, you may need to put `mpir.dll` and `mpfr.dll` into the main project folder (where the C++ files are)
-   Select 'Release' and 'x64' in the toolbar and build using Ctrl+Shift+B

## How to run?
Currently there is no way to run at custom parameters other than compiling yourself. But just to prove that this actually works, I uploaded a compiled binary of the current state of the code as a release here on GitHub.
In its current state the program will render the image above.