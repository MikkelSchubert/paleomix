import nox

nox.options.sessions = [
    "style",
    "tests",
]


@nox.session
def style(session: nox.Session) -> None:
    session.install("black~=23.9.1")
    session.run("black", "--check", ".")


@nox.session
def lints(session: nox.Session) -> None:
    session.install("ruff==0.0.292")
    session.run(
        "ruff",
        "paleomix",
        "tests",
    )


@nox.session()
def typing(session: nox.Session) -> None:
    session.install(".")
    session.install("pytest~=7.4")
    session.install("pyright==1.1.329")
    session.run(
        "pyright",
        "paleomix",
        "tests",
    )


@nox.session()
def tests(session: nox.Session) -> None:
    session.install(".")
    session.install("pytest~=7.4")
    session.install("coverage[toml]~=7.3")
    session.install("pytest-cov~=4.1")

    session.run(
        "pytest",
        "tests",
        # TODO: Treat warnings as errors
        # "-Werror",
        "--cov",
        "paleomix",
        "--cov",
        "tests",
        "--cov-report=xml",
        "--cov-report=term-missing",
        "-m",
        "not slow",
        "--quiet",
    )


@nox.session(python=["3.7", "3.8", "3.9", "3.10", "3.11"])
def full_tests(session: nox.Session) -> None:
    session.install(".")
    session.install("pytest~=7.4")

    session.run("pytest", "tests", "--quiet")
