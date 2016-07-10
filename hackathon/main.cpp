#include <time.h>
#include <stdlib.h>
#include <stdio.h>
#include <stdarg.h>
#include <fstream>
#include <string.h>
#include "utils.h"
#include <FTGL/ftgl.h>

#define pi 3.141592

#define SCREEN_W 2560
#define SCREEN_H 1440
#define PLOT_W 0.3*1.35
#define PLOT_H 0.3*1.20*SCREEN_W/SCREEN_H
#define NB_TICKS 5


#define W_TO_H_LENGTH(W) (W*(float)SCREEN_W/(float)SCREEN_H)
#define H_TO_W_LENGTH(H) (H*(float)SCREEN_H/(float)SCREEN_W)



#define GL_TO_FT_COORD_X(X) ( 0.5*(X+1.0)*SCREEN_W )
#define GL_TO_FT_COORD_Y(Y) ( 0.5*(Y+1.0)*SCREEN_H )

// The coordinate system used for the window is:
//          +y 
//          |
//          |
//          |
// -x<------------->+x
//          |
//          |
//          |
//          -y
// In other words the top left corner of the window as coordinates (x,y) = (-1.0,+1.0)
// while the bottom right corner as coordinates (x,y) = (+1.0,-1.0).

double L = 1.;              	// heat exchanger length (meters)
double dx = 0.1;           	// grid size for spatial integration (meters)
int Nx = L/dx;			      	// number of length steps


double rho_i = 1000.;       	// density of inner fluid (kg/m3)
double rho_o = 1000.;       	// density of outer fluid (kg/m3)

double cp_i = 4200.;        	// specific heat capacity of inner fluid (j/kg.K)
double cp_o = 4200.;        	// specific heat capacity of outer fluid (j/kg.K)

double Q_i = 1e-4;          	// flowrate of inner fluid (m3/s)
double Q_o = 1e-4;          	// flowrate of outer fluid (m3/s)

double r_i = 1e-1;          	// radius of inner pipe (meters)
double r_o = 2e-1;          	// radius of outer pipe (meters)

double A_i = pi * pow(r_i, 2);
double A_o = pi * (pow(r_o, 2) - pow(r_i, 2));

double u_i = Q_i / A_i;
double u_o = Q_o / A_o;

double h = 1000.;          		// heat transfer coefficient of heat exchanger (j/m2.K.s)

double a_i = Q_i * rho_i * cp_i / (A_i * rho_i);
double a_o = Q_o * rho_o * cp_o / (A_o * rho_o);

double b_i = 2 * pi * r_i * h / (A_i * rho_i);
double b_o = 2 * pi * r_i * h / (A_o * rho_i);

double dt = 5e-5;

double T0_i = 50.;
double T0_o = 10.;

void integrate(double * Ti, double * To, int Nx){
	double * Ti_copy = (double*)malloc(sizeof(double) * Nx);
	double * To_copy = (double*)malloc(sizeof(double) * Nx);
	memcpy(Ti_copy, Ti, sizeof(double) * Nx); 
	memcpy(To_copy, To, sizeof(double) * Nx);
	
	for (int j = 1; j < Nx; j++){
        Ti[j] = Ti_copy[j] + dt * (-a_i/dx * (Ti_copy[j] - Ti_copy[j-1])
        - b_i * (Ti_copy[j] - To_copy[j]));       
       
        To[j] = To_copy[j] + dt * (-a_o/dx * (To_copy[j] - To_copy[j-1])
        + b_o * (Ti_copy[j] - To_copy[j]));
	}

	free(Ti_copy);
	free(To_copy);
}

void integrate2(double * Ti, double * To, int Nx){
	double * Ti_copy = (double*)malloc(sizeof(double) * Nx);
	double * To_copy = (double*)malloc(sizeof(double) * Nx);
	memcpy(Ti_copy, Ti, sizeof(double) * Nx); 
	memcpy(To_copy, To, sizeof(double) * Nx);
	
	for (int j = 1; j < Nx; j++){
        Ti[j] = Ti_copy[j] + dt * (-a_i/dx * (Ti_copy[j] - Ti_copy[j-1])
        - b_i * (Ti_copy[j] - To_copy[j]));
	}
       	for (int j = 0; j < Nx-1; j++){
        To[j] = To_copy[j] + dt * (+a_o/dx * (To_copy[j+1] - To_copy[j])
        + b_o * (Ti_copy[j] - To_copy[j]));
	}

	free(Ti_copy);
	free(To_copy);
}

SDL_Surface *convertSurface(SDL_Surface *surface)
{
  Uint32 rmask = 0, gmask = 0, bmask = 0, amask = 0;
  /*
    Change byte order (depends if the system uses Big endian or not)
    necessery for multiple platform programming.
  */
#if SDL_BYTEORDER == SDL_BIG_ENDIAN

  rmask = 0xff000000;
  gmask = 0x00ff0000;
  bmask = 0x0000ff00;
  amask = 0x000000ff;
#else

  rmask = 0x000000ff;
  gmask = 0x0000ff00;
  bmask = 0x00ff0000;
  amask = 0xff000000;
#endif

  SDL_PixelFormat format = *(surface->format);
  format.BitsPerPixel = 32;
  format.BytesPerPixel = 4;
  format.Rmask = rmask;
  format.Gmask = gmask;
  format.Bmask = bmask;
  format.Amask = amask;

  return SDL_ConvertSurface(surface, &format, SDL_SWSURFACE);
}

GLuint loadImage(const char *filename)
{
  GLuint id;
  SDL_Surface *surface = IMG_Load(filename);
  SDL_Surface *surface32 = convertSurface(surface);

  glEnable(GL_TEXTURE_2D);
  glEnable(GL_BLEND);
  glGenTextures(1, &id);
  glBindTexture(GL_TEXTURE_2D, id);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
  glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
  SDL_LockSurface(surface32);
  glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA, surface32->w, surface32->h, 0, GL_RGBA, GL_UNSIGNED_BYTE, surface32->pixels);
  SDL_UnlockSurface(surface32);
  glBindTexture(GL_TEXTURE_2D, 0);
  SDL_FreeSurface(surface);
  SDL_FreeSurface(surface32);
  return id;
}


void drawImage(GLuint id, float x, float y, float w, float h)
{
  glLoadIdentity();
  glColor4f(1.0, 10, 1.0, 1.0);
  glTranslated(x, y, 0.0);
  glBindTexture(GL_TEXTURE_2D, id);
  glBegin(GL_QUADS);
  glTexCoord2f(0.0, 1.0);
  glVertex2f(-0.5*w, -0.5*h);
  glTexCoord2f(1.0, 1.0);
  glVertex2f(+0.5*w, -0.5*h);
  glTexCoord2f(1.0, 0.0);
  glVertex2f(+0.5*w, 0.5*h);
  glTexCoord2f(0.0, 0.0);
  glVertex2f(-0.5*w, 0.5*h);
  glEnd();
  glBindTexture(GL_TEXTURE_2D, 0);
}

void drawArray(double *array, int nbElems, float min, float max, float r, float g, float b)
{
  float tempMin = FLT_MAX, tempMax = 0.0; 
  float x = -0.5*PLOT_W;

  //  glLoadIdentity();
  glColor3f(r, g, b);
  glLineWidth(4.0);
  glBegin(GL_LINES);
  for(int i = 0 ; i < nbElems-1 ; i++)
    {
      glVertex2f(x, PLOT_H/(max-min)*array[i] + 0.5*PLOT_H-PLOT_H*max/(max-min));
      glVertex2f(x+PLOT_W / (float)(nbElems-1), PLOT_H/(max-min)*array[i+1] + 0.5*PLOT_H-PLOT_H*max/(max-min));
      x += PLOT_W / (float)(nbElems-1);
    }
  glEnd();
}

void drawPlot(double *Ti, double *To, int Nx, FTGLPixmapFont *font, float posX, float posY, float minX, float maxX, float minY, float maxY)
{
  static GLuint yAxis = loadImage("y-axis.png");
  
  glLoadIdentity();
  glTranslated(posX, posY, 0.0);
  glColor3f(0.0, 0.0, 0.0);
  glLineWidth(2.0);
  glBegin(GL_LINES);
  glVertex2f(-0.5*PLOT_W, -0.5*PLOT_H);
  glVertex2f(+0.5*PLOT_W, -0.5*PLOT_H);
  
  glVertex2f(-0.5*PLOT_W, -0.5*PLOT_H);
  glVertex2f(-0.5*PLOT_W, +0.5*PLOT_H);

  glVertex2f(-0.5*PLOT_W, +0.5*PLOT_H);
  glVertex2f(+0.5*PLOT_W, +0.5*PLOT_H);

  glVertex2f(+0.5*PLOT_W, -0.5*PLOT_H);
  glVertex2f(+0.5*PLOT_W, +0.5*PLOT_H);
  glEnd();

  glBegin(GL_LINES);
  for(int i = 0 ; i < NB_TICKS ; i++)
    {
      glColor3f(0.0, 0.0, 0.0);
      glVertex2f(-0.5*PLOT_W-0.01, -0.5*PLOT_H+(i+1)*PLOT_H/(float)(NB_TICKS+1));
      glVertex2f(-0.5*PLOT_W-0.0, -0.5*PLOT_H+(i+1)*PLOT_H/(float)(NB_TICKS+1));

      glVertex2f(-0.5*PLOT_W+(i+1)*PLOT_W/(float)(NB_TICKS+1), -0.5*PLOT_H-0.01);
      glVertex2f(-0.5*PLOT_W+(i+1)*PLOT_W/(float)(NB_TICKS+1), -0.5*PLOT_H-0.0);   
    }
  glEnd();
  
  glEnable(GL_LINE_STIPPLE);
  glLineStipple(1, 0x0F);
  glBegin(GL_LINES);
  for(int i = 0 ; i < NB_TICKS ; i++)
    {
      glColor3f(0.8, 0.8, 0.8);
      glVertex2f(-0.5*PLOT_W, -0.5*PLOT_H+(i+1)*PLOT_H/(float)(NB_TICKS+1));
      glVertex2f(+0.5*PLOT_W, -0.5*PLOT_H+(i+1)*PLOT_H/(float)(NB_TICKS+1));

      glVertex2f(-0.5*PLOT_W+(i+1)*PLOT_W/(float)(NB_TICKS+1), -0.5*PLOT_H);
      glVertex2f(-0.5*PLOT_W+(i+1)*PLOT_W/(float)(NB_TICKS+1), +0.5*PLOT_H);      
    }
  glEnd();
  glDisable(GL_LINE_STIPPLE);
  
  glPixelTransferf(GL_RED_BIAS, -1.0f);
  glPixelTransferf(GL_GREEN_BIAS, -1.0f);
  glPixelTransferf(GL_BLUE_BIAS, -1.0f);

  char tick[10];
  for(int i = 0 ; i < NB_TICKS + 2 ; i++)
    {
      sprintf(tick, "%.2f", minX + i*(maxX-minX)/(NB_TICKS+1));
      (*font).Render(tick, -1, FTPoint(GL_TO_FT_COORD_X(posX - 0.550*PLOT_W + i*PLOT_W/(NB_TICKS+1)), GL_TO_FT_COORD_Y(posY - 0.5*PLOT_H-0.05), 0.0));
    }
  for(int i = 0 ; i < NB_TICKS + 2 ; i++)
    {
      sprintf(tick, "%.1f", minY + i*(maxY-minY)/(NB_TICKS+1));
      (*font).Render(tick, -1, FTPoint(GL_TO_FT_COORD_X(posX - 0.5*PLOT_W-0.05), GL_TO_FT_COORD_Y(posY - 0.510*PLOT_H + i*PLOT_H/(NB_TICKS+1)), 0.0));
    }

  (*font).Render("Length of the Pipe (m)", -1, FTPoint(GL_TO_FT_COORD_X(posX-0.09), GL_TO_FT_COORD_Y(posY-0.65*PLOT_H), 0.0));
  drawImage(yAxis, posX - 0.5*PLOT_W-0.09, posY+0.0, 57.0*2.0/SCREEN_W, 240.0*2.0/SCREEN_H);

  glLoadIdentity();
  glTranslated(posX, posY, 0.0);  
  drawArray(Ti, Nx, T0_o, T0_i, 1.0, 0.0, 0.0);
  drawArray(To, Nx, T0_o, T0_i, 0.0, 0.0, 1.0);
}

int main()
{
  bool quit = false, rerun = false, startSimulation = false, isLengthChange = false;
  SDL_Window *window;
  SDL_Event event;
  double * Ti = (double*)malloc(sizeof(double) * Nx);
  double * To = (double*)malloc(sizeof(double) * Nx);
  double * Ti2 = (double*)malloc(sizeof(double) * Nx);
  double * To2 = (double*)malloc(sizeof(double) * Nx);
  float gui_rho_i = rho_i, gui_rho_o = rho_o, gui_cp_i = cp_i, gui_cp_o = cp_o, gui_q_i = Q_i, gui_q_o = Q_o, gui_t0_i = T0_i, gui_t0_o = T0_o, gui_l = L;
  int nbIntegration = 0;
  for (int j = 0; j < Nx; j++){
    Ti[j] = T0_i;
    To[j] = T0_o;
    Ti2[j] = T0_i;
    To2[j] = T0_o;
  }
  
  srand(time(NULL));
  
  window = initWindow("HEAT", SCREEN_W, SCREEN_H);
  FTGLPixmapFont font("OpenSans-Bold.ttf");
  font.FaceSize(24);

  GLuint imageID = loadImage("DoubleExchanger_Transparent.png");
  GLuint imageID2 = loadImage("DoubleExchanger_Transparent2.png");
  while (!quit)
    {
      while(SDL_PollEvent(&event))
        {
  	  processEvent(&event);
  	  if(event.type == SDL_QUIT)
            {
  	      quit = true;
            }
  	  if(event.type == SDL_KEYDOWN)
  	    {
  	      switch(event.key.keysym.sym)
  		{
  		case SDLK_ESCAPE:
  		  quit = true;
  		  break;
  		}
  	    }
        }
      windowGUINewFrame(window);
      glClearColor(1.0f, 1.0f, 1.0f, 1.f);
      glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
      //>>>>>> GUI
      ImGui::Begin("Parameters", NULL, ImGuiWindowFlags_AlwaysAutoResize);
      ImGui::SliderFloat("tube density (kg/m^3)", &gui_rho_i, 500, 1000, "%4.0f");
      //      ImGui::SameLine(0.0, 104);
      ImGui::SliderFloat("shell density (kg/m^3)", &gui_rho_o, 500, 1000, "%4.0f");
      ImGui::SliderFloat("tube heat capacity (J/kg.C)", &gui_cp_i, 500, 5000, "%4.0f");
      //      ImGui::SameLine();
      ImGui::SliderFloat("shell heat capacity (J/kg.C)", &gui_cp_o, 500, 5000, "%4.0f");
      ImGui::SliderFloat("tube flow rate (m^3/s)", &gui_q_i, 1e-6, 1e-3, "%.7f");
      //      ImGui::SameLine(0.0, 88);
      ImGui::SliderFloat("shell flow rate (m^3/s)", &gui_q_o, 1e-6, 1e-3, "%.7f");
      ImGui::SliderFloat("tube inlet temperature (C)", &gui_t0_i, 0, 100, "%3.0f");
      //      ImGui::SameLine(0.0, 24);
      ImGui::SliderFloat("shell inlet temperature (C)", &gui_t0_o, 0, 100, "%3.0f");
      if(isLengthChange == false)
	{
	  isLengthChange = ImGui::SliderFloat("tube length (m)", &gui_l, 0.1, 2, "%3.1f");
	}
      else
	{
	  ImGui::SliderFloat("tube length (m)", &gui_l, 0.1, 2, "%3.1f");
	}
      rerun = ImGui::Button("simulate");
      ImGui::End();
      ImGui::Begin("Time");
      ImGui::Text("%f (s)", (float)(nbIntegration*dt));
      ImGui::End();
      ImGui::Render();
      //<<<<<< GUI
      drawImage(imageID, -0.6, -0.55, 0.75, 0.5);
      drawImage(imageID2, +0.55, -0.6, 0.8, 0.95);
      drawPlot(Ti, To, Nx, &font, -0.69, +0.35, 0, L, T0_o, T0_i);
      drawPlot(Ti2, To2, Nx, &font, +0.69, +0.35, 0, L, T0_o, T0_i);
      SDL_GL_SwapWindow(window);

      if(rerun)
	{
	  nbIntegration = 0;
	  startSimulation = true;

	  if(isLengthChange)
	    {
	      L = gui_l;
	      Nx = L/dx;
	      free(Ti);
	      free(To);
	      Ti = (double*)malloc(sizeof(double) * Nx);
	      To = (double*)malloc(sizeof(double) * Nx);
	      free(Ti2);
	      free(To2);
	      Ti2 = (double*)malloc(sizeof(double) * Nx);
	      To2 = (double*)malloc(sizeof(double) * Nx);
	      isLengthChange = false;
	    }
	  
	  T0_i = gui_t0_i;
	  T0_o = gui_t0_o;
	  for (int j = 0; j < Nx; j++)
	    {
	      Ti[j] = T0_i;
	      To[j] = T0_o;
	      Ti2[j] = T0_i;
	      To2[j] = T0_o;
	    }

	  rho_i = gui_rho_i;
	  rho_o = gui_rho_o;
	  cp_i = gui_cp_i;
	  cp_i = gui_cp_o;
	  Q_i = gui_q_i;
	  Q_o = gui_q_o;

	  A_i = pi * pow(r_i, 2);
	  A_o = pi * (pow(r_o, 2) - pow(r_i, 2));
	  a_i = Q_i * rho_i * cp_i / (A_i * rho_i);
	  a_o = Q_o * rho_o * cp_o / (A_o * rho_o);
	  b_i = 2 * pi * r_i * h / (A_i * rho_i);
	  b_o = 2 * pi * r_i * h / (A_o * rho_i);
	}
      if(startSimulation)
	{
	  for(int i = 0 ; i < 10 ; i++)
	    {
	      nbIntegration++;
	      integrate(Ti, To, Nx);
	      integrate2(Ti2, To2, Nx);
	    }
	}
    }

  free(Ti);
  free(To);
  free(Ti2);
  free(To2);
  glDeleteTextures(1, &imageID);
  glDeleteTextures(1, &imageID2);
  freeWindow(window);
  return 0;
}
