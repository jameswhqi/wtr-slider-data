# Data and code for the paper "An accurate and efficient measure of welfare tradeoff ratios"

## Data and analyses

Anonymized raw data files are `data/expt{1,2,3}.json` for the three experiments, respectively.

### Installing dependencies

The following has been tested on macOS. It should be similar on Linux, while Windows users might need to make some adjustments.

To run the data analyses, first set up an R environment. We highly recommend using [Nix](https://nixos.org/) to install the dependencies. After [installing Nix](https://nixos.org/download), run `nix develop` in this repository. You might need to add
```
experimental-features = nix-command flakes
```
to your [`nix.conf`](https://nixos.org/manual/nix/stable/command-ref/conf-file) file. After a while, you will be dropped in a development shell with all the dependencies installed and on `$PATH`. If you are using VS Code, you can instead run the command `nix develop -c sh -c 'code -n .'` in a terminal (not the integrated terminal in VS Code) to open a window with the environment variables properly set.

If you don't want to use Nix, you can also install the dependencies manually. See [`flake.nix`](flake.nix) for the required R packages.

### Running analyses

We use [`targets`](https://docs.ropensci.org/targets/) to manage the pipelines of data processing, model building, and output generation. To build all the analyses, `cd` to the `data/` directory and run `./build.R`. This will take some time and computation, mostly due to sampling in Stan. The built artifacts will be written to `store-expt{1,2,3}/`, which you can examine using methods described below. The output data used for producing the figures in the paper are created by targets with the suffix `Output` and will be written to `output/`, which you can mostly ignore.

To build or examine the results of individual targets, create an R terminal and run `library(targets)` and then `Sys.setenv(TAR_PROJECT = "expt1")`, which specifies that you are dealing with expt1 (similarly for expt2 and expt3). Then run `tar_manifest()` to see the list of available targets for that experiment. Run `tar_make(<target_name>)` to build an individual target and all its dependencies. Run `tar_read(<target_name>)` to read (as a variable) the result of an individual target (it must have been built).

`data/script-expt{1,2,3}.R` define the list of targets in each experiment. `tar_target(<target_name>, ...)` creates individual targets, while other top-level functions such as `tar_stanFit()` are defined in `data/R/common.R` and create multiple targets at once.

For any model fit with raw RStan (defined by `tar_stanFit()`), there is usually a target with the suffix `Fit` containing the fit object and a target with the suffix `Summary` containing the summary of relevant parameters. For any model fit with `brms` (usually defined by `tar_target()` itself), the target has the suffix `Fit` containing the fit object, which, when printed, by default shows the posterior means instead of medians. To see the medians, use `print(<fit_object>, robust = T)`.

## Experiment interfaces

`experiments/expt{1,2,3}/` contain the web pages of the three experiments, which can also be found at `https://experiments.evullab.org/qi-games-{2,4,7}/` (note the set notation; pick one number out of the three).

If the links no longer work or if you want to test data writing, you can run an experiment locally by `cd`ing into one of the subdirectories (e.g., `expt1/`) and running `python3 -m http.server --cgi` (assuming a relatively new version of Python is installed). Then the web page can be viewed at `http://localhost:8000/`. Data will be written to the `data/` subdirectory at the end of the experiment.

## Standalone Lambda Slider

`standalone/` contains a standalone version of the Lambda Slider, which can also be found at `https://experiments.evullab.org/lambda-slider/`.

There are 6 URL parameters you can add:

| Parameter | Meaning | Value type | Default value |
| - | - | - | - |
| `selfName` | Name of self | string | `You` |
| `otherName` | Name of target | string | `Other` |
| `selfText` | Text after self name | string | `%20receive` |
| `otherText` | Text after target name | string | `%20receives` |
| `init` | Initial slider position in terms of λ | number | (Random within [-2, 2]) |
| `selfOnTop` | Whether self payoff is on top | `true` or `false` | `true` |

Example: `https://experiments.evullab.org/lambda-slider/?otherName=Alice&init=0.5&selfOnTop=false`

If the link no longer works, you can run `python3 -m http.server` and view it at `http://localhost:8000/`.

### Embedding the standalone Lambda Slider in Qualtrics

In Qualtrics, create a Text / Graphic question (Content type: Text), and edit the HTML of the text to be something like:
```html
<div>This is the description of the question.</div>
<div>
  <iframe src="https://experiments.evullab.org/lambda-slider/?otherName=Alice&init=0.5" width="100%" height="260"></iframe>
</div>
```
where "Alice" and "0.5" can be replaced by piped text like `${e://Field/someField}` if necessary.

Whenever the slider is moved, the `iframe` will call (in JavaScript) something like `window.parent.postMessage({lambda: 1.5})` to send the current λ to the parent window, so you can edit the JavaScript of the question to be something like:
```javascript
let listener;

Qualtrics.SurveyEngine.addOnload(function() {
  this.disableNextButton();
  listener = e => {
    if (e.data.lambda) {
      Qualtrics.SurveyEngine.setEmbeddedData("lambda", e.data.lambda);
      this.enableNextButton();
    }
  };
  window.addEventListener("message", listener);
});

Qualtrics.SurveyEngine.addOnUnload(function() {
  window.removeEventListener("message", listener);
});
```
such that the participant has to move the slider at least once before they can go to the next page, and the response will be recorded in the `lambda` field of the embedded data (you can change it to other fields).

## Questions

Please submit an issue if you have any questions or have difficulty reproducing the results.
