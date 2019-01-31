.. _testing:

Testing
################

EXOSIMS provides both unit tests for component level input/output and functionality validation as well as an end-to-end test suite for full scale operational tests.

End to End Testing
=====================

End to end testing is provided by the ``e2eTests.py`` method in the main EXOSIMS directory. This utility is intended to be executed directly (i.e. ``python e2eTests.py``) and looks for test scripts in the ``Scripts/TestScripts/`` directory in the EXOSIMS folder hierarchy.

For each test script, ``e2eTests`` will:

* Create a ``MissionSim`` object
* Execute a simulation via ``MissionSim.run_sim()``
* Reset the simulation via ``MissionSim.reset_sim()``
* Execute a second simulation from the same object via ``MissionSim.run_sim()``

The test suite will record a ``PASS`` or ``FAIL`` condition for each individual step, and will print a summary of results for all scripts at the end of execution. 

There are currently 14 test scripts as of 01/23/2019
* 01_all_defaults.json - Fundamental prototype module loading and execution test
* 02_KnownRV_FAP=1_WFIRSTObs_staticEphem.json - Testing KnownRV with static Ephemeris
* 03_EarthTwin_Coronagraph_GarrettComp.json - Testing Earth Twin Planet Population and Garrett Completeness
* 04_KeplerLike_Occulter_linearJScheduler.json - Testing KeplerLike, Occulter, and LinearJScheduler
* 05_SAG13_Forecaster_BrownCompleteness.json - Testing BrownCompleteness, SAG13, and Forecaster
* 06_Occulter.json - Testing Occulter in a quasi relistic universe
* 07_TKtestingOBfromFile0.json - Testing loading OB from file SampleOB.csv with OBendTime[-1]=380days and 200days total, termination by exceeding MissionLife=1yr and MissonPortion=1
* 08_TKtestingOBfromFile1.json - Testing termination by exceeding exoplanetObsTime with MissionLife=1yr and MissonPortion=0.5479452
* 09_TKtestingOBfromFile2.json - Testing termination by exceeding OBendTime[-1] using SampleOB1.csv
* 10_TKtestingOBfromFile3.json - Testing termination by exceeding MissionLife using SampleOB2.csv
* 11_TKtestingOBfromFile4.json - Testing termination by exceeding exoplanetObsTime with MissionLife=1yr and MissionPortion=0.5
* 12_TKtestingOBfromJSON.json - Testing generation of Observing Blocks (OB) from JSON script, OBduration=14days, MissionPortion=0.5, MissionLife=1yr
* 13_TKtestingSingleOBFullPortion.json - Testing OB are automatically generated for MissionLife=1yr and MissionPortion=1
* 14_TKtestingSingleOBHalfPortion.json - Testing OB are automatically generated for MissionLife=1yr and MissionPortion=0.5


Unit Testing
====================

Unit tests are implemented using Python's ``unittest`` framework (see https://docs.python.org/2/library/unittest.html). The unit tests are in the ``tests`` directory under the EXOSIMS root directory.  The ``tests`` directory contains the same folder hierarchy as the EXOSIMS directory, with a separate folder for each module type, a folder for module prototypes, a folder for utility method unit tests and a ``TestSupport`` folder for test-specific code and input scripts. There are three types of unit tests:

* In folder ``ModuleName`` the ``test_ModuleName.py`` object is used to test basic functionality of all implementations of that module (including the prototype).  Test coverage in these objects is limited only to methods found the prototype, and is only executed on specifically overloaded methods in module implementations.
* In folder ``ModuleName``, any ``test_ImplementationName.py`` object is used to test specific functionality of _only_ that module implementation.  These tests are only needed to test specific expected behavior of implementations, or methods in implementations that are not overloading a prototype method.
* In folder ``Prototypes``, the ``test_ModuleName.py`` object is used to test specific functionality of prototype implementations. 

Unit tests can be executed by running:
:: 
    python -m unittest discover -v

from the EXOSIMS root directory.  This will execute all available unit tests (but not the end to end tests). Individual tests can be executed by running:
::
    python -m unittest -v tests.ModuleName.testName

to run all tests within that object or:
::
    python -m unittest -v tests.ModuleName.testName.testMethodName
    
to run a single individual test.

