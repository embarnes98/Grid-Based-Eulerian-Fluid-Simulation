#include <SFML/Graphics.hpp> 
#include "Fluid-Sim-SFML.h"
#include <cstdlib>
#include <iostream>

int main()
{
  // Set up fluid grid instance
  auto diffusion = 1e-6;
  auto viscosity = 0;
  std::unique_ptr<FluidGrid> fluidGrid(new FluidGrid(diffusion, viscosity));

  // SFML window/drawing stuff
  auto windowScale = 5;
  sf::RenderWindow window(sf::VideoMode(windowScale*Nx, windowScale*Ny), "SFML works!", sf::Style::Close | sf::Style::Titlebar);
  window.setView(sf::View(sf::FloatRect(0, 0, Nx, Ny)));
  auto numColorValues = 4 * Nx * Ny;
  auto *pixels = new sf::Uint8[numColorValues];
  int i, j, pixelNum, colorValueNum = 0;

  while (colorValueNum < numColorValues)
  {
    pixels[colorValueNum++] = 0; // R
    pixels[colorValueNum++] = 0; // G
    pixels[colorValueNum++] = 0; // B
    pixels[colorValueNum++] = 255; // A
  }

  sf::Image image;
  sf::Texture texture;
  sf::Sprite sprite;

  // Text and graphics stuff
  sf::Font font;
  font.loadFromFile("arial.ttf");
  auto fontSize = 10;
  sf::Text frameRateTxt, leftEdgeTxt, rightEdgeTxt, topEdgeTxt, bottomEdgeTxt;

  frameRateTxt.setFont(font);
  frameRateTxt.setCharacterSize(fontSize);
  frameRateTxt.setFillColor(sf::Color::Yellow);

  auto leftEdge = fluidGrid->getLeftEdge();
  auto rightEdge = fluidGrid->getRightEdge();
  auto topEdge = fluidGrid->getTopEdge();
  auto bottomEdge = fluidGrid->getBottomEdge();

  bool viewOptimisationGraphics = true;
  
  // Time stuff
  
  sf::Clock clock;
  float t0, t1, t_updateGraphics, dt;
  t0 = clock.getElapsedTime().asSeconds();
  t_updateGraphics = t0;
  bool timeToUpdateGraphics;
  auto speedFactor = 10;

  // Mouse stuff
  sf::Vector2i mousePos(Nx/(2*windowScale), Ny/(2*windowScale));
  sf::Vector2i prevMousePos;
  int dvx, dvy, jmin, jmax, imin, imax;
  auto blotRadius = 5;
  auto velocityFactor = 1;

  // Simulation loop

  while(window.isOpen())
  {
    prevMousePos = mousePos;
    mousePos = sf::Mouse::getPosition(window)/windowScale;
    mousePos.x = std::max(blotRadius + 1, std::min(mousePos.x, Nx - blotRadius - 2));
    mousePos.y = std::max(blotRadius + 1, std::min(mousePos.y, Ny - blotRadius - 2));
    sf::Event event;
    while (window.pollEvent(event))
    {
      if (event.type == sf::Event::Closed)
        window.close();
      if (sf::Keyboard::isKeyPressed(sf::Keyboard::Escape))
        window.close();
      if (sf::Mouse::isButtonPressed(sf::Mouse::Left))
      {
        dvx = (mousePos.x - prevMousePos.x)*velocityFactor;
        dvy = (mousePos.y - prevMousePos.y)*velocityFactor;
        for (j = mousePos.y - blotRadius; j <= mousePos.y + blotRadius; j++)
        {
          for (i = mousePos.x - blotRadius; i <= mousePos.x + blotRadius; i++)
          {
            if (IN_CIRCLE(i, j, mousePos.x, mousePos.y, blotRadius))
            {
              fluidGrid->addDensity(i, j, 100);
              fluidGrid->addVelocity(i, j, dvx, dvy);
              fluidGrid->updateLocalGrid(i, j); 
            }
          }
        }
      }
    }

    t1 = clock.getElapsedTime().asSeconds();
    timeToUpdateGraphics = false; 
    if (t1 - t_updateGraphics > .33f)
    {
      t_updateGraphics = t1;
      timeToUpdateGraphics = true;
    }
    dt = t1 - t0;
    t0 = t1;
    if (timeToUpdateGraphics)
    {
      frameRateTxt.setString(std::to_string(int(1/dt)) + "fps");
    }
    dt *= speedFactor;
    fluidGrid->step(dt);

    if (fluidGrid->isOptimised() && viewOptimisationGraphics)
    {
      for (colorValueNum = 0, pixelNum = 0; colorValueNum < numColorValues; colorValueNum += 4, pixelNum++)
      {
        pixels[colorValueNum    ] = std::max(30 * (1 - fluidGrid->getSurfaceFlag(pixelNum)), fluidGrid->getDensity(pixelNum));
        pixels[colorValueNum + 1] = fluidGrid->getDensity(pixelNum);
        pixels[colorValueNum + 2] = std::max(30 * fluidGrid->getSurfaceFlag(pixelNum), fluidGrid->getDensity(pixelNum));
      }
      if (timeToUpdateGraphics)
      {
        leftEdge = fluidGrid->getLeftEdge();
        rightEdge = fluidGrid->getRightEdge();
        topEdge = fluidGrid->getTopEdge();
        bottomEdge = fluidGrid->getBottomEdge();
      }
      if (leftEdge != rightEdge && topEdge != bottomEdge)
      {
        for (i = leftEdge; i <= rightEdge; i++)
        {
          pixels[IX(i, topEdge   ) * 4 + 1] = 50;
          pixels[IX(i, bottomEdge) * 4 + 1] = 50;
        }
        for (i = topEdge; i <= bottomEdge; i++)
        {
          pixels[IX(leftEdge, i ) * 4 + 1] = 50;
          pixels[IX(rightEdge, i) * 4 + 1] = 50;
        }
      }
    }
    else
    {
      for (colorValueNum = 0, pixelNum = 0; colorValueNum < numColorValues; colorValueNum += 4, pixelNum++)
      {
        pixels[colorValueNum    ] = fluidGrid->getDensity(pixelNum);
        pixels[colorValueNum + 1] = fluidGrid->getDensity(pixelNum);
        pixels[colorValueNum + 2] = fluidGrid->getDensity(pixelNum);
      }
    }
    
    window.clear();
    image.create(Nx, Ny, pixels);
    texture.loadFromImage(image);
    sprite.setTexture(texture, true);
    window.draw(sprite);
    window.draw(frameRateTxt);
    window.display();
  }

  delete [] pixels;

  return EXIT_SUCCESS;
}