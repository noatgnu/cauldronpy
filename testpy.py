import click

@click.command()
@click.option("--boole", "-b", help="Boolean option", default=False)
def main(boole: bool):
    print(boole)

if __name__ == "__main__":
    main()