#!/usr/bin/env python
# -*- coding: utf-8 -*-

from setuptools import setup

requirements = [
    'numpy >= 1.12.1',
    'scipy >= 0.19.0',
    'click >= 8.0',
]

test_requirements = [
    'pytest'
]

setup(
    name = 'superfocus',
    version = '1.8',
    description = 'SUPER-FOCUS: A tool for agile functional analysis of shotgun metagenomic data',
    author = 'Genivaldo G.Z. Silva',
    author_email = 'genivaldo.gueiros@gmail.com',
    url = 'https://github.com/metageni/super-focus',
    packages = [
        'superfocus_app',
    ],
    package_dir={'superfocus_app': 'superfocus_app'},
    include_package_data=True,
    install_requires=requirements,
    zip_safe=False,
    keywords='superfocus_app',
    classifiers=[
        'Development Status :: 2 - Pre-Alpha',
        'Intended Audience :: Developers',
        'Natural Language :: English',
        'Programming Language :: Python :: 3.6',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
    ],
    test_suite='tests',
    tests_require=test_requirements,
    entry_points={
        'console_scripts': [
            'superfocus = superfocus_app.cli:main',
        ]
    },
)
