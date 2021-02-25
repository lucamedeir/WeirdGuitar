#include<iostream>
#include <vector>
#include<gsl/gsl_sf_trig.h>
#include<SDL2/SDL.h>

using namespace std::string_view_literals;

//Screen dimension constants
const int SCREEN_WIDTH = 640;
const int SCREEN_HEIGHT = 480;

//The window we'll be rendering to
SDL_Window* gWindow = NULL;

//The window renderer
SDL_Renderer* gRenderer = NULL;

//Current displayed texture
SDL_Texture* gTexture = NULL;

bool init()
{
    //Initialization flag
    bool success = true;

    //Initialize SDL
    if (SDL_Init(SDL_INIT_VIDEO) < 0)
    {
        printf("SDL could not initialize! SDL_Error: %s\n", SDL_GetError());
        success = false;
    }
    else
    {
        //Create window
        gWindow = SDL_CreateWindow("Weird Guitar", SDL_WINDOWPOS_UNDEFINED, SDL_WINDOWPOS_UNDEFINED, SCREEN_WIDTH, SCREEN_HEIGHT, SDL_WINDOW_SHOWN);
        if (gWindow == NULL)
        {
            printf("Window could not be created! SDL Error: %s\n", SDL_GetError());
            success = false;
        }
        else
        {
            //Create renderer for window
            gRenderer = SDL_CreateRenderer(gWindow, -1, SDL_RENDERER_ACCELERATED);
            if (gRenderer == NULL)
            {
                printf("Renderer could not be created! SDL Error: %s\n", SDL_GetError());
                success = false;
            }
            else
            {
                //Initialize renderer color
                SDL_SetRenderDrawColor(gRenderer, 0xFF, 0xFF, 0xFF, 0xFF);

            }
        }
    }

    return success;
}


void close()
{
    //Free loaded image
    SDL_DestroyTexture(gTexture);
    gTexture = NULL;

    //Destroy window    
    SDL_DestroyRenderer(gRenderer);
    SDL_DestroyWindow(gWindow);
    gWindow = NULL;
    gRenderer = NULL;

    //Quit SDL subsystems
    SDL_Quit();
}

int main(int argc, char *argv[])
{
    //Main loop flag
    bool quit = false;

    //Event handler
    SDL_Event e;

    std::vector<double> X;
#define X_SIZE 100
    X.reserve(X_SIZE);
    const double dx = 1.0 / (double)X_SIZE;

    //populate vector X
    for (int i = 0; i < X_SIZE; ++i) {
        X.push_back(i*dx);
    }

    std::vector<double> Y;
    Y.reserve(X_SIZE);

    const double k = 10.0;

    for (auto& x : X)
    {
        Y.push_back(gsl_sf_sin(k*x));
    }

    //Start up SDL and create window
    if (!init())
    {
        printf("Failed to initialize!\n");
    }
    else
    {
        const double omega = 10.0;
        const double dt = 0.0005;
        double t = 0.0;
        long count = 0;

        //While application is running
        while (!quit)
        {
            //Handle events on queue
            while (SDL_PollEvent(&e) != 0)
            {
                //User requests quit
                if (e.type == SDL_QUIT)
                {
                    quit = true;
                }
            }
            //Clear screen
            SDL_SetRenderDrawColor(gRenderer, 0xFF, 0xFF, 0xFF, 0xFF);
            SDL_RenderClear(gRenderer);


            //Draw blue horizontal line
            SDL_SetRenderDrawColor(gRenderer, 0x00, 0x00, 0xFF, 0xFF);

            t += dt;

            auto xi = X.begin();
            auto yi = Y.begin();

            while (xi != X.end() - 1) {
                auto x1 = *xi;
                auto x2 = *(xi+1);
                auto y1 = *yi;
                auto y2 = *(yi + 1);

                *yi = gsl_sf_sin(k*x1 + t * omega);

                SDL_RenderDrawLine(gRenderer,
                    (int)SCREEN_WIDTH*x1,
                    SCREEN_HEIGHT / 2 - (int)SCREEN_HEIGHT / 2 * y1,
                    (int)SCREEN_WIDTH * x2,
                    SCREEN_HEIGHT / 2 - (int)SCREEN_HEIGHT / 2 * y2);

                

                xi++;
                yi++;
            }

            

            //Update screen
            SDL_RenderPresent(gRenderer);

            count++;
        }
    }

    //Free resources and close SDL
    close();

	return 0;
}