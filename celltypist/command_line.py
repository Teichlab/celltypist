import os
import math
import click
from celltypist import logger, annotate, defaults, models


def show_banner():
    logger.info(r"""
                    oooo  oooo      .                           o8o               .
                    `888  `888    .o8                           `"'             .o8
 .ooooo.   .ooooo.   888   888  .o888oo oooo    ooo oo.ooooo.  oooo   .oooo.o .o888oo
d88' `"Y8 d88' `88b  888   888    888    `88.  .8'   888' `88b `888  d88(  "8   888
888       888ooo888  888   888    888     `88..8'    888   888  888  `"Y88b.    888
888   .o8 888    .o  888   888    888 .    `888'     888   888  888  o.  )88b   888 .
`Y8bod8P' `Y8bod8P' o888o o888o   "888"     .8'      888bod8P' o888o 8""888P'   "888"
                                        .o..P'       888
                                        `Y8P'       o888o""")


def show_config(config: dict):
    logger.info(f"Configuration:")
    logger.info(f"\t-Input: {config['indata']}")
    logger.info(f"\t-Model: {config['model']}")
    logger.info(f"\t-Output prefix: {config['outpref'] if config['outpref'] else '(none)'}")
    logger.info(f"\t-Output path: {config['outdir']}")
    logger.info(f"\t-Cell count: {config['cell_count']}")
    logger.info(f"\t-Chunk count: {config['chunk_count']} (each one of {config['chunk_size']} cells)")
    logger.info(f"\t-CPUs: {config['cpus']}")


def show_help_and_exit(message: str):
    ctx = click.get_current_context()
    click.echo(click.style(message, fg="red"))
    click.echo()
    ctx.fail(ctx.get_help())


def write_xlsx(result, outpref, outdir):
    output_filename = f"{outpref}annotation_result.xlsx"
    result.write_excel(
        os.path.join(outdir, output_filename))


def write_all_csv_files(result, outpref, outdir):
    # write summary
    summary_filename = f"{outpref}summary.csv"
    result.summary.to_csv(
        os.path.join(outdir, summary_filename))
    # write labels
    labels_filename = f"{outpref}predicted_labels.csv"
    result.predicted_labels.to_csv(
        os.path.join(outdir, labels_filename))
    # write prob matrix
    probability_filename = f"{outpref}probability_matrix.csv"
    result.probability_matrix.to_csv(
        os.path.join(outdir, probability_filename))


@click.command()
@click.option("--indata", help="Input csv matrix (cells by genes). Gene IDs should be Gene Names.", type=str)
@click.option("--model", default=defaults.model, help="Model used to make the predictions.", type=str)
@click.option("--outpref", default="", help="Output prefix for all output files.", type=str)
@click.option("--outdir", default="", help="Output directory for all output files.", type=str)
@click.option("--xlsx", is_flag=True, default=False, help="Merge output files into a single XLSX.")
@click.option("--chunk", default=defaults.chunk_size, help="Chunk sizes to read (adjust for memory performance).", type=int)
@click.option("--cpus", default=defaults.max_cpus, help="Limit the numbre of CPUs (default uses all avaiable).", type=int)
@click.option("--quiet", is_flag=True, default=False, help="Hide all console output.")
def main(indata: str, model: str, outpref: str, outdir: str, xlsx: bool, chunk: int, cpus: int, quiet: bool):

    # validate input file
    if not indata or not os.path.exists(indata):
        show_help_and_exit(f"Missing or invalid input file: '{indata}'")

    # validate model name/file
    if model not in models.get_all_models() and not os.path.exists(model):
        show_help_and_exit(f"Missing or invalid model: '{model}'. Avaiable models are: {','.join(models.get_all_models())}")

    if not outdir:
        outdir = os.getcwd()
        logger.debug(f"No output directory provided. Using current directory: {os.getcwd()}")

    with open(indata) as fh:
        total_size = sum(1 for line in fh)

    config = {
        "indata": indata,
        "model": model,
        "outpref": outpref,
        "outdir": outdir,
        "chunk_size": chunk,
        "chunk_count": math.ceil(total_size/chunk),
        "cell_count": total_size,
        "cpus": cpus,
        "quiet": quiet
    }

    if not quiet:
        show_banner()
    if not quiet:
        show_config(config)

    result = annotate(
        filename=config["indata"],
        model=config["model"],
        chunk_size=config["chunk_size"],
        cpus=config["cpus"],
        quiet=config["quiet"])

    if xlsx:
        write_xlsx(result, config["outpref"], config["outdir"])
    else:
        write_all_csv_files(result, config["outpref"], config["outdir"])
