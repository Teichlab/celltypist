import os
import click
from . import logger, defaults, models
from .annotate import annotate

def show_banner():
    logger.info(r"""
                    oooo  oooo      .                           o8o               .
                     888   888    .o8                           `"'             .o8
 .ooooo.   .ooooo.   888   888  .o888oo oooo    ooo oo.ooooo.  oooo   .oooo.o .o888oo
d88'  "Y8 d88   88b  888   888    888     88.  .8    888   88b  888  d88(  "8   888
888       888ooo888  888   888    888      88..8     888   888  888  `"Y88b.    888
888   .o8 888    .o  888   888    888 .    `888'     888   888  888  o.  )88b   888 .
 Y8bod8P  `Y8bod8P' o888o o888o   "888"     .8'      888bod8P' o888o 8""888P'   "888"
                                        .o..P'       888
                                        `Y8P'       o888o             @Teichmann lab""")


def show_config(config: dict):
    logger.info(f"🛠️ Configuration:")
    for key, value in config.items():
        logger.info(f"\t📥 {key}: {value}")


def show_help_and_exit(message: str):
    ctx = click.get_current_context()
    click.echo(click.style(message, fg="red"))
    click.echo()
    ctx.fail(ctx.get_help())


def write_xlsx(result, prefix, outdir):
    output_filename = f"{prefix}annotation_result.xlsx"
    result.write_excel(
        os.path.join(outdir, output_filename))


def write_all_csv_files(result, prefix, outdir):
    #labels
    labels_filename = f"{prefix}predicted_labels.csv"
    result.predicted_labels.to_csv(
        os.path.join(outdir, labels_filename))
    #prob table
    probability_filename = f"{prefix}probability_matrix.csv"
    result.probability_table.to_csv(
        os.path.join(outdir, probability_filename))

@click.command()
@click.option("-i", "--indata", help="Input count matrix (e.g., .csv) or Scanpy object (.h5ad). Genes should be provided as gene symbols.", type=str)
@click.option("-m", "--model", default="", help="Model used for predictions. Default to using the 'Immune_All_Low.pkl' model.", type=str)
@click.option("--transpose-input", is_flag=True, default=False, help="Transpose the input matrix if --indata file is provided in the gene-by-cell format. Note Celltypist needs a cell-by-gene matrix as input.")
@click.option("--majority-voting", is_flag=True, default=False, help="Run the majority voting classifier to refine the predicted labels.")
@click.option("-oc", "--over-clustering", default='auto', help="Input file with clustering result of one cell per line, or a string key specifying an existing metadata column in the AnnData. Default to using a heuristic over-clustering approach based on input data size. Ignored if --majority-voting is not set.", type=str)
@click.option("-o", "--outdir", default="", help="Directory to store the output file/files. Default to the current working directory.", type=str)
@click.option("--xlsx", is_flag=True, default=False, help="Merge output files into a single Excel (.xlsx).")
@click.option("-p", "--prefix", default="", help="Prefix for the output file/files. Default to no specific prefix.", type=str)
@click.option("--update-models", is_flag=True, default=False, help="Download the latest models from the remote server.")
@click.option("--quiet", is_flag=True, default=False, help="Hide the banner and configure information during the run.")
def main(indata: str, model: str, transpose_input: bool, majority_voting: bool, over_clustering,
         outdir: str, xlsx: bool, prefix: str, update_models: bool, quiet: bool):

    #update models or not
    if update_models:
        models.update_models()
        exit(0)

    #validate input file
    if not indata or not os.path.exists(indata):
        show_help_and_exit(f"🛑 Missing or invalid input file: '{indata}'")

    #validate model
    if not model:
        model = models.get_default_model()
        logger.info(f"🔖 No model provided. Using the deault: '{model}'")
    if model not in models.get_all_models() and not os.path.exists(model):
        show_help_and_exit(f"🛑 Invalid model name: '{model}'. Available models are: {', '.join(models.get_all_models())}")

    #output dir
    if not outdir:
        outdir = os.getcwd()
        logger.warn(f"👀 No output directory provided. Using the current directory: {outdir}")

    #config settings
    if not majority_voting:
        config = {
                "indata": indata,
                "model": model,
                "transpose-input": transpose_input,
                "majority-voting": majority_voting,
                "outdir": outdir,
                "xlsx": xlsx,
                "prefix": prefix,
                "quiet": quiet
                }
    else:
        config = {
                "indata": indata,
                "model": model,
                "transpose-input": transpose_input,
                "majority-voting": majority_voting,
                "over-clustering": over_clustering,
                "outdir": outdir,
                "xlsx": xlsx,
                "prefix": prefix,
                "quiet": quiet
                }

    #quiet or not
    if not quiet:
        show_banner()
        show_config(config)

    #celltyping and majority voting
    if over_clustering == 'auto':
        over_clustering = None
    result = annotate(
        filename=indata,
        model=model,
        transpose_input=transpose_input,
        majority_voting=majority_voting,
        over_clustering=over_clustering)

    #write output
    if xlsx:
        write_xlsx(result, config["prefix"], config["outdir"])
    else:
        write_all_csv_files(result, config["prefix"], config["outdir"])
