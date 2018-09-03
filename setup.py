
from setuptools import setup, find_packages


def main():
    setup(
        name="excitingparser",
        version="0.1",
        description="NOMAD parser implementation for Exciting.",
        license="APACHE 2.0",
        package_dir={'': 'parser/parser-exciting'},
        packages=find_packages(),
        install_requires=[
            'nomadcore'
        ],
    )

if __name__ == "__main__":
    main()
