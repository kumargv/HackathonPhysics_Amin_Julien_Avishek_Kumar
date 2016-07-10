#!/bin/sh
g++ -std=c++11 `sdl2-config --cflags` `pkg-config --cflags --libs-only-L --libs-only-l ftgl` -I imgui $C_INCLUDE_PATH main.cpp imgui/imgui*.cpp `sdl2-config --libs` -lSDL2_image -lGL -lGLEW -o main
