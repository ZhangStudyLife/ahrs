AHRS: Attitude and Heading Reference Systems
============================================

Welcome! ``ahrs`` is an open source Python toolbox for attitude estimation
using the most known algorithms, methods and resources.

It is designed to be flexible and simple to use, making it a great option for
fast prototyping, testing and integration with your own Python projects.

This package collects functions and utilities to help you understand and use
the most common techniques for attitude estimation, and in no way it is
recommended to be used commercially.

All algorithms and implementations have their proper documentation and
references, in case you need further clarification of their usage.

New in version 0.4
------------------

- New attitude estimation algorithms: FKF and UKF.
- New submodules: ``sensors`` and ``geodesy``.
- More methods for Quaternion and Direction Cosine Matrix representations.
- Improved documentation and examples.
- All functions and objects are type hinted.

For a complete list of changes, please refer to the `changelog
<https://github.com/Mayitzin/ahrs/blob/master/CHANGELOG.md>`_.

Attitude Estimation
-------------------

The available algorithms to compute attitude from inertial / magnetic sensors
are:

=============  =========  =============  ============
Algorithm      Gyroscope  Accelerometer  Magnetometer
=============  =========  =============  ============
AQUA           Optional   YES            Optional
Complementary  YES        YES            Optional
Davenport's    NO         YES            YES
EKF            YES        YES            YES
FAMC           NO         YES            YES
FLAE           NO         YES            YES
FKF            YES        YES            YES
Fourati        YES        YES            YES
FQA            NO         YES            Optional
Integration    YES        NO             NO
Madgwick       YES        YES            Optional
Mahony         YES        YES            Optional
OLEQ           NO         YES            YES
QUEST          NO         YES            YES
ROLEQ          NO         YES            YES
SAAM           NO         YES            YES
Tilt           NO         YES            Optional
TRIAD          NO         YES            YES
UKF            YES        YES            Optional
=============  =========  =============  ============

Deprecations
------------

Submodules ``io`` and ``plot`` are removed. Loading and visualizing the data is
left to the preference of the user.

.. toctree::
   :maxdepth: 1
   :hidden:

   installation
   geodesy
   attitude_representations
   filters
   metrics
   sensors
   constants
   nomenclature
   bibliography

Indices
=======

* :ref:`genindex`
* :ref:`search`
