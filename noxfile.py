import nox

nox.options.default_venv_backend = "uv|virtualenv"
nox.options.sessions = [
    "style",
    "lints",
    "tests",
]


SOURCES = (
    "noxfile.py",
    "src",
    "tests",
    "typings",
)


@nox.session
def style(session: nox.Session) -> None:
    session.install("--group", "linting")

    # Replaces `black --check`
    session.run("ruff", "format", "--check", *SOURCES)
    # Replaces `isort --check-only`
    session.run("ruff", "check", "--select", "I", *SOURCES)


@nox.session
def lints(session: nox.Session) -> None:
    session.install("--group", "linting")

    session.run("ruff", "check", *SOURCES)


@nox.session()
def typing(session: nox.Session) -> None:
    session.install("-e", ".", "--group", "typing")

    session.run("basedpyright", *SOURCES)


@nox.session()
def tests(session: nox.Session) -> None:
    # Install in development mode to for coverage analysis
    session.install("-e", ".", "--group", "testing")

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


@nox.session(python=["3.10", "3.11", "3.12", "3.13", "3.14"])
def full_tests(session: nox.Session) -> None:
    session.install("-e", ".", "--group", "testing")

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
