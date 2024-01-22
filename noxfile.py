import nox

nox.options.sessions = [
    "style",
    "tests",
]


SOURCES = (
    "noxfile.py",
    "paleomix",
    "setup.py",
    "tests",
    "typings",
)


class Requirements:
    COVERAGE = "coverage[toml]~=7.3"
    NOX = "nox~=2023.4.22"
    PYRIGHT = "pyright==1.1.338"
    PYTEST = "pytest~=7.4"
    PYTEST_COV = "pytest-cov~=4.1"
    RUFF = "ruff==0.1.14"


@nox.session
def style(session: nox.Session) -> None:
    session.install(Requirements.RUFF)
    # Replaces `black --check`
    session.run("ruff", "format", "--check", *SOURCES)
    # Replaces `isort --check-only`
    session.run("ruff", "check", "--select", "I", *SOURCES)


@nox.session
def lints(session: nox.Session) -> None:
    session.install(Requirements.RUFF)
    session.run("ruff", "check", *SOURCES)


@nox.session()
def typing(session: nox.Session) -> None:
    session.install(".")
    session.install(Requirements.PYTEST)
    session.install(Requirements.NOX)
    session.install(Requirements.PYRIGHT)
    session.run("pyright", *SOURCES)


@nox.session()
def tests(session: nox.Session) -> None:
    session.install("-e", ".")
    session.install(Requirements.PYTEST)
    session.install(Requirements.COVERAGE)
    session.install(Requirements.PYTEST_COV)

    session.run(
        "python3",
        # Run tests in development mode (enables extra checks)
        "-X",
        "dev",
        "-m",
        "pytest",
        "tests",
        "--cov",
        "paleomix",
        "--cov",
        "tests",
        "--cov-report=xml",
        "--cov-report=term-missing",
        "--no-cov-on-fail",
        # Exclude slow tests, including tests that run the pipeline via subprocess
        "-m",
        "not slow",
        "--quiet",
        # Treat warnings (deprections, etc.) as errors
        "-Werror",
        # Re-run failed tests, or all tests if there were no failures
        "--last-failed",
        "--last-failed-no-failures",
        "all",
    )


@nox.session(python=["3.7", "3.8", "3.9", "3.10", "3.11", "3.12"])
def full_tests(session: nox.Session) -> None:
    session.install(".")
    session.install("pytest~=7.4")

    session.run(
        "python3",
        # Run tests in development mode (enables extra checks)
        "-X",
        "dev",
        "-m",
        "pytest",
        "tests",
        "--quiet",
        # Treat warnings (deprections, etc.) as errors
        "-Werror",
    )
