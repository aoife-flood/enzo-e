
***********************
Particle Merging
***********************

Merging Star Particles Tests
############################

These tests run an idealised setup of a uniform spherical distribution of star particles, which collapse under their own gravity and merge together. They are run using either gravity solvers or a constant radial infall velocity applied to all particles. The purpose of these tests is to check if momentum and mass are conserved throughout the merging process.

Parameters
**********

Initial

 * Infall_speed: Constant radial velocity directed towards the centre of collapse that is applied to all particles

 * centre: The centre of collapse and the centre of the sphere.

 * drift_velocity: Additional velocity that can be added to the star particles.

 * truncation_radius: The radius of the sphere of star paticles.

Method 

 * list = [ "pm_deposit", "gravity", "pm_update","merge_stars"];

   * To turn off gravity remove "gravity" and "pm_deposit" from list.

   * To turn off merging remove "merge_stars" from list.


 * merging_radius_cells: Determines how close star particles have to get to merge, the numerical value refers to cells rather than any physical unit.


How to run tests
****************

* Change the parameters file.
  
* Run the script run.sh with "nohup ./run.sh".
  
  This will:
  
  * Remove DIR directories and .png files from previous runs
    
  * Run Enzo-e for the specified parameters
    
  * Run make_images.py in parallel

  * Run data.py

  * Run multiple.py
    
* After running it will output:
  
  * Images of a cross section of the star particles at each time interval. (Output of make_images.py)
    
    .. figure:: cross_section_image.png
          :width: 590px
          :align: center
          :height: 470px
	  :alt: alternate text
	  :figclass: align-center

	  Cross Section of Sphere of star paricles at initial time
		     
  * A .txt file named merge_data_1_0.txt (where 1_0 refers to a merging cells radius of 1.0), which contains the sum of the velcities and the mass of all the star particles as well as the number of particles for each time interval. (output of data.py)

    The data in the file is in the form:
    
    +--------------+--------------+--------------+------+------------------+------+
    | X velocities | Y velocities | Z velocities | Mass | No. of particles | Time |
    +--------------+--------------+--------------+------+------------------+------+
    
  * A graph of Mass, Normalised momentum and No. of particles VS Time/Collapse Time. (output of data.py)
    Momentum in these graphs is the magnitude of the total momentum normalised by dividing by (G*M^3)/R or M*S depending on whether gravity or infall speed is used, where G is the universal gravitational constant, M is the total mass, R is the radius and s is infall speed.

    .. figure:: Mass_momentum_particles_graph_0_3centreofblock.png
          :width: 590px
          :align: center
          :height: 470px
          :alt: alternate text
          :figclass: align-center
		     
          Graph of Mass, Momentum and No. of Particles VS Time for 0.3 Merging radius with Gravity Off
  * A graph of Mass, Normalised momentum and No. of particles VS Time/Collapse Time for different merge radii using data from previous runs



Test 1
******

For the first test, merging is turned off to check that momentum and mass are conserved.

* Remove "merge_stars" from Method: list in the parameters file.

* Check that "pm _deposit" and "gravity" are in Method: list.

* Check that Initial: Infall_speed=0.0
  
* Run test using "nohup ./run.sh"

.. figure:: Graph1.png
    :width: 590px
    :align: center
    :height: 470px
    :alt: alternate text
    :figclass: align-center

    Graph of Momentum VS Time with merging off

For the rest of the tests merging is turned on, so "merge_stars" should be added back into Method: list in the parameters file.

Test 2
******

For this test gravity is turned off, the centre of the collapse is positioned in the centre of a block to ensure any errors are not coming from errors in particles being copied across blocks, and the truncation radius is made very small so that there are fewer particles. Momentum and mass should be conserved, the particle number should decrease.

* Remove "pm _deposit" and "gravity" from Method: list in the parameters file.

* Set Initial: infall_speed to a suitable value that will allow all particles to reach the centre by the end of the stopping time.

* Set collapse_centre in Initial to be [3.086e24, 3.0856e24, 3.086e24] in parameters file.

* Set upper/lower bounds in Domain to be [12.344e24,12.344e24,12.344e24] and [-12.344e24,-12.344e24,-12.344e24] in parameters file.

* Set truncation_radius in Initial to be 3.0e23 in parameters file.

* Run test using "nohup ./run.sh"

* Run the test for multiple merge radii by changing merging_cell_radius in Method: merge_stars

* Graph all radii on one plot by running multiple.py

.. figure:: Test2.png
    :width: 590px
    :align: center
    :height: 470px
    :alt: alternate text
    :figclass: align-center

    Graph of Momentum, Mass and No. of particles VS Time with gravity off and small truncation radius in one block

Test 3
******

Test 3 is like test 2 but with a larger truncation radius and more particles. It should show similar results to test 2. Momentum and mass should be conserved, the particle number should decrease.

* Set truncation_radius in Initial to be 3.086e24 in parameters file.

* Run test using "nohup ./run.sh"

* Run the test for multiple merge radii by changing merging_cell_radius in Method: merge_stars

* Graph all radii on one plot by running multiple.py

.. figure:: Test3.png
    :width: 590px
    :align: center
    :height: 470px
    :alt: alternate text
    :figclass: align-center

    Graph of Momentum, Mass and No. of particles VS Time with gravity off and large truncation radius in one block

Test 4
******

For this test, the same set up is used, but the collapse centre is changed so that the collapse and merging will take place across blocks. If the results of this test differ greatly from the previous test it will mean there is a problem occuring when particles are being copied across blocks. Momentum and mass should be conserved, the particle number should decrease.

* Set collapse_centre in Initial to be [0.0,0.0,0.0] in parameters file.

* Set upper/lower bounds in Domain to be [6.172e24, 6.172e24, 6.172e24] in parameters file.

* Run test using "nohup ./run.sh"

* Run the test for multiple merge radii by changing merging_cell_radius in Method: merge_stars

* Graph all radii on one plot by running multiple.py

.. figure:: Test4.png
    :width: 590px
    :align: center
    :height: 470px
    :alt: alternate text
    :figclass: align-center

    Graph of Momentum, Mass and No. of particles VS Time with gravity off and large truncation radius


Test 5
******

For this test gravity is turned back on, momentum and mass should still be conserved and the results should be similar to the previous test.

* Add "pm _deposit" and "gravity" from Method: list in the parameters file.

* Run test using "nohup ./run.sh"

* Run the test for multiple merge radii by changing merging_cell_radius in Method: merge_stars

* Graph all radii on one plot by running multiple.py

.. figure:: Test5.png
    :width: 590px
    :align: center
    :height: 470px
    :alt: alternate text
    :figclass: align-center

    Graph of Momentum, Mass and No. of particles VS Time with gravity on and large truncation radius


Test 6
******

For this the drift velocity is changed to a non-zero number to check that momentum and mass are still conserved.

* In the parameters file set Initial: drift_velocity = 

* Run test using "nohup ./run.sh"

* Run the test for multiple merge radii by changing merging_cell_radius in Method: merge_stars

* Graph all radii on one plot by running multiple.py

.. figure:: Test6.png
    :width: 590px
    :align: center
    :height: 470px
    :alt: alternate text
    :figclass: align-center

    Graph of Momentum, Mass and No. of particles VS Time with gravity on, large truncation radius, and non-zero drift velocity

