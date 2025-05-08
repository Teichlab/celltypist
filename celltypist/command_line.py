import os
import click
from . import logger, models
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
    logger.info(f"⚙️ Configuration:")
    for key, value in config.items():
        logger.info(f"\t🛠️ {key}: {value}")


def show_help_and_exit(message: str):
    ctx = click.get_current_context()
    click.echo(click.style(message, fg="red"))
    click.echo()
    ctx.fail(ctx.get_help())


@click.command()
@click.option("-i", "--indata", help="Path to the input count matrix (.csv/txt/tsv/tab/mtx) or AnnData (.h5ad). Genes should be provided as gene symbols.", type=click.Path(exists=True, dir_okay=False))
@click.option("-m", "--model", default=None, help="Model used for predictions. If not provided, default to using the `Immune_All_Low.pkl` model.", type=str)
@click.option("--transpose-input", is_flag=True, default=False, help="Transpose the input matrix if `-i / --indata` file is provided in the gene-by-cell format. Note Celltypist requires the cell-by-gene format.")
@click.option("-gf", "--gene-file", default=None, type=click.Path(exists=False), help="Path to the file which stores each gene per line corresponding to the genes used in the provided mtx file. Ignored if `-i / --indata` is not provided in the mtx format.")
@click.option("-cf", "--cell-file", default=None, type=click.Path(exists=False), help="Path to the file which stores each cell per line corresponding to the cells used in the provided mtx file. Ignored if `-i / --indata` is not provided in the mtx format.")
@click.option("-mo", "--mode", default="best_match", help="Choose the cell type with the largest score/probability as the final prediction (`best_match`), or enable a multi-label classification (`prob_match`), which assigns 0 (i.e., unassigned), 1, or >=2 cell type labels to each query cell.", type=click.Choice(['best_match','prob_match']), show_default=True)
@click.option("-ps", "--p-thres", default=0.5, help="Probability threshold for the multi-label classification. Ignored if `--mode` is `best_match`.", type=float, show_default=True)
@click.option("--majority-voting", is_flag=True, default=False, help="Refine the predicted labels by running the majority voting classifier after over-clustering.")
@click.option("-oc", "--over-clustering", default='auto', help="Input file with the over-clustering result of one cell per line, or a string key specifying an existing metadata column in the AnnData object. If not provided, default to using a heuristic over-clustering approach according to the size of input data. Ignored if `--majority-voting` is not set.", type=str, show_default=True)
@click.option("--use-GPU", is_flag=True, default=False, help="Whether to use GPU for over clustering on the basis of `rapids-singlecell`. Ignored if `--majority-voting` is not set.")
@click.option("-mp", "--min-prop", default=0, help="For the dominant cell type within a subcluster, the minimum proportion of cells required to support naming of the subcluster by this cell type. Ignored if `--majority-voting` is not set.", type=float, show_default=True)
@click.option("-o", "--outdir", default=None, help="Directory to store the output files and (if `--plot-results` is set) figures. Default to the current working directory.", type=click.Path(exists=False))
@click.option("-p", "--prefix", default="", help="Prefix for the output files and (if `--plot-results` is set) figures. Default to no prefix used.", type=str)
@click.option("--xlsx", is_flag=True, default=False, help="Merge output tables into a single Excel (.xlsx).")
@click.option("--plot-results", is_flag=True, default=False, help="Plot the prediction results as multiple figures as well.")
@click.option("--update-models", is_flag=True, default=False, help="Download the latest models from the remote server.")
@click.option("--show-models", is_flag=True, default=False, help="Show all the available models and their descriptions.")
@click.option("--quiet", is_flag=True, default=False, help="Hide the banner and configuration information during the run.")
def main(indata: str, model: str, transpose_input: bool, gene_file: str, cell_file: str, mode: str, p_thres: float, majority_voting: bool, over_clustering: str, use_GPU: bool, min_prop: float,
         outdir: str, prefix: str, xlsx: bool, plot_results: bool, update_models: bool, show_models: bool, quiet: bool):
    """Celltypist: a tool for semi-automatic cell type annotation"""

    #update models or not
    if update_models:
        models.download_models(force_update=True)
        exit(0)

    #show all models
    if show_models:
        md = models.models_description()
        for _, row in md.iterrows():
            row = row.tolist()
            logger.info(row[0] + '   ' + row[1])
        exit(0)

    #validate model
    if model is None:
        model = models.get_default_model()
        logger.info(f"🔖 No model provided. Using the deault: '{model}'")
    if not os.path.isfile(model) and model not in models.get_all_models():
        show_help_and_exit(f"🛑 Invalid model name: '{model}'. Available models are: {', '.join(models.get_all_models())}")

    #output dir
    if outdir is None:
        outdir = os.getcwd()
        logger.warn(f"👀 No output directory provided. Using the current directory: '{outdir}'")
    if not os.path.isdir(outdir):
        show_help_and_exit(f"🛑 Output directory '{outdir}' does not exist")

    #config settings
    config = {
                "indata": indata,
                "model": model,
                "transpose-input": transpose_input,
                "gene-file": gene_file,
                "cell-file": cell_file,
                "mode": mode,
                "p-thres": p_thres,
                "majority-voting": majority_voting,
                "outdir": outdir,
                "prefix": prefix,
                "xlsx": xlsx,
                "plot-results": plot_results,
                "quiet": quiet
             }
    if majority_voting:
        config["over-clustering"] = over_clustering
        config["use-GPU"] = use_GPU
        config["min-prop"] = min_prop

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
        gene_file=gene_file,
        cell_file=cell_file,
        mode = mode.replace('_', ' '),
        p_thres = p_thres,
        majority_voting=majority_voting,
        over_clustering=over_clustering,
        use_GPU = use_GPU,
        min_prop = min_prop)

    #write output
    result.to_table(folder = outdir, prefix = prefix, xlsx = xlsx)

    #plot result
    if plot_results:
        result.to_plots(folder = outdir, prefix = prefix, plot_probability = True)

if __name__ == "__main__":
    main()
