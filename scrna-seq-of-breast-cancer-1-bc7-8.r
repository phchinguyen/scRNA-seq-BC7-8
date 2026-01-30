{
 "cells": [
  {
   "cell_type": "code",
   "execution_count": 1,
   "id": "82c668e5",
   "metadata": {
    "_execution_state": "idle",
    "_uuid": "051d70d956493feee0c6d64651c6a088724dca2a",
    "execution": {
     "iopub.execute_input": "2026-01-30T04:26:10.325585Z",
     "iopub.status.busy": "2026-01-30T04:26:10.323107Z",
     "iopub.status.idle": "2026-01-30T04:26:11.328325Z",
     "shell.execute_reply": "2026-01-30T04:26:11.326375Z"
    },
    "papermill": {
     "duration": 1.014375,
     "end_time": "2026-01-30T04:26:11.331588",
     "exception": false,
     "start_time": "2026-01-30T04:26:10.317213",
     "status": "completed"
    },
    "tags": []
   },
   "outputs": [],
   "source": [
    "# Download 10x Genomics Dataset - \n",
    "download.file(\n",
    "  \"https://cf.10xgenomics.com/samples/cell-exp/9.0.0/320k_scFFPE_16-plex_GEM-X_FLEX_BreastCancer1_BC7-8/320k_scFFPE_16-plex_GEM-X_FLEX_BreastCancer1_BC7-8_count_sample_raw_feature_bc_matrix.h5\",\n",
    "  destfile = \"320k_scFFPE_BreastCancer1_BC7_raw_feature_bc_matrix.h5\",\n",
    "  mode = \"wb\"\n",
    ")"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 2,
   "id": "ba1c5a2e",
   "metadata": {
    "execution": {
     "iopub.execute_input": "2026-01-30T04:26:11.371214Z",
     "iopub.status.busy": "2026-01-30T04:26:11.338822Z",
     "iopub.status.idle": "2026-01-30T04:26:45.928581Z",
     "shell.execute_reply": "2026-01-30T04:26:45.924770Z"
    },
    "papermill": {
     "duration": 34.597717,
     "end_time": "2026-01-30T04:26:45.932296",
     "exception": false,
     "start_time": "2026-01-30T04:26:11.334579",
     "status": "completed"
    },
    "tags": []
   },
   "outputs": [
    {
     "name": "stderr",
     "output_type": "stream",
     "text": [
      "Loading required package: SeuratObject\n",
      "\n"
     ]
    },
    {
     "name": "stderr",
     "output_type": "stream",
     "text": [
      "Loading required package: sp\n",
      "\n"
     ]
    },
    {
     "name": "stderr",
     "output_type": "stream",
     "text": [
      "\n",
      "Attaching package: ‘SeuratObject’\n",
      "\n",
      "\n"
     ]
    },
    {
     "name": "stderr",
     "output_type": "stream",
     "text": [
      "The following objects are masked from ‘package:base’:\n",
      "\n",
      "    intersect, t\n",
      "\n",
      "\n"
     ]
    },
    {
     "name": "stderr",
     "output_type": "stream",
     "text": [
      "\n",
      "Attaching package: ‘dplyr’\n",
      "\n",
      "\n"
     ]
    },
    {
     "name": "stderr",
     "output_type": "stream",
     "text": [
      "The following objects are masked from ‘package:stats’:\n",
      "\n",
      "    filter, lag\n",
      "\n",
      "\n"
     ]
    },
    {
     "name": "stderr",
     "output_type": "stream",
     "text": [
      "The following objects are masked from ‘package:base’:\n",
      "\n",
      "    intersect, setdiff, setequal, union\n",
      "\n",
      "\n"
     ]
    },
    {
     "name": "stderr",
     "output_type": "stream",
     "text": [
      "Warning message:\n",
      "“Feature names cannot have underscores ('_'), replacing with dashes ('-')”\n"
     ]
    },
    {
     "name": "stderr",
     "output_type": "stream",
     "text": [
      "Warning message:\n",
      "“Feature names cannot have underscores ('_'), replacing with dashes ('-')”\n"
     ]
    }
   ],
   "source": [
    "library(Seurat)\n",
    "library(dplyr)\n",
    "library(patchwork)\n",
    "# Read and Create Seurat object\n",
    "sc_data <- Read10X_h5(\"320k_scFFPE_BreastCancer1_BC7_raw_feature_bc_matrix.h5\")\n",
    "sc_obj <- CreateSeuratObject(counts = sc_data, project = \"320k_scFFPE\")\n",
    "\n",
    "sc <- CreateSeuratObject(counts= sc_data, project = \"320k_scFFPE\", min.cells = 3, min.features = 200)"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 3,
   "id": "b67d9056",
   "metadata": {
    "execution": {
     "iopub.execute_input": "2026-01-30T04:26:45.943157Z",
     "iopub.status.busy": "2026-01-30T04:26:45.941068Z",
     "iopub.status.idle": "2026-01-30T04:26:45.966435Z",
     "shell.execute_reply": "2026-01-30T04:26:45.964713Z"
    },
    "papermill": {
     "duration": 0.032958,
     "end_time": "2026-01-30T04:26:45.968568",
     "exception": false,
     "start_time": "2026-01-30T04:26:45.935610",
     "status": "completed"
    },
    "tags": []
   },
   "outputs": [
    {
     "data": {
      "text/plain": [
       "An object of class Seurat \n",
       "18811 features across 34171 samples within 1 assay \n",
       "Active assay: RNA (18811 features, 0 variable features)\n",
       " 1 layer present: counts"
      ]
     },
     "metadata": {},
     "output_type": "display_data"
    }
   ],
   "source": [
    "sc"
   ]
  },
  {
   "cell_type": "markdown",
   "id": "a5acf10e",
   "metadata": {
    "papermill": {
     "duration": 0.003088,
     "end_time": "2026-01-30T04:26:45.974826",
     "exception": false,
     "start_time": "2026-01-30T04:26:45.971738",
     "status": "completed"
    },
    "tags": []
   },
   "source": [
    "# Standard Preprocessing\n",
    "## 1. Quality Control"
   ]
  },
  {
   "cell_type": "markdown",
   "id": "0a526172",
   "metadata": {
    "papermill": {
     "duration": 0.003619,
     "end_time": "2026-01-30T04:26:45.982110",
     "exception": false,
     "start_time": "2026-01-30T04:26:45.978491",
     "status": "completed"
    },
    "tags": []
   },
   "source": [
    "## 2. Normalization"
   ]
  },
  {
   "cell_type": "markdown",
   "id": "ff61c1cc",
   "metadata": {
    "papermill": {
     "duration": 0.003036,
     "end_time": "2026-01-30T04:26:45.988326",
     "exception": false,
     "start_time": "2026-01-30T04:26:45.985290",
     "status": "completed"
    },
    "tags": []
   },
   "source": [
    "## 3. Important Features"
   ]
  },
  {
   "cell_type": "markdown",
   "id": "67d49cc9",
   "metadata": {
    "papermill": {
     "duration": 0.002992,
     "end_time": "2026-01-30T04:26:45.994277",
     "exception": false,
     "start_time": "2026-01-30T04:26:45.991285",
     "status": "completed"
    },
    "tags": []
   },
   "source": [
    "## 4. PCA"
   ]
  },
  {
   "cell_type": "markdown",
   "id": "742fd743",
   "metadata": {
    "papermill": {
     "duration": 0.003023,
     "end_time": "2026-01-30T04:26:46.000231",
     "exception": false,
     "start_time": "2026-01-30T04:26:45.997208",
     "status": "completed"
    },
    "tags": []
   },
   "source": [
    "## 5. Dimensionality of Dataset"
   ]
  },
  {
   "cell_type": "markdown",
   "id": "cf5f02f9",
   "metadata": {
    "papermill": {
     "duration": 0.002914,
     "end_time": "2026-01-30T04:26:46.006049",
     "exception": false,
     "start_time": "2026-01-30T04:26:46.003135",
     "status": "completed"
    },
    "tags": []
   },
   "source": [
    "### a. Heatmap"
   ]
  },
  {
   "cell_type": "markdown",
   "id": "21d314fe",
   "metadata": {
    "papermill": {
     "duration": 0.002998,
     "end_time": "2026-01-30T04:26:46.011931",
     "exception": false,
     "start_time": "2026-01-30T04:26:46.008933",
     "status": "completed"
    },
    "tags": []
   },
   "source": [
    "### b. Elbow plot method"
   ]
  },
  {
   "cell_type": "markdown",
   "id": "0ce06866",
   "metadata": {
    "papermill": {
     "duration": 0.004346,
     "end_time": "2026-01-30T04:26:46.021889",
     "exception": false,
     "start_time": "2026-01-30T04:26:46.017543",
     "status": "completed"
    },
    "tags": []
   },
   "source": [
    "## 6. Cell Clustering"
   ]
  },
  {
   "cell_type": "markdown",
   "id": "c057bd2e",
   "metadata": {
    "papermill": {
     "duration": 0.003051,
     "end_time": "2026-01-30T04:26:46.028675",
     "exception": false,
     "start_time": "2026-01-30T04:26:46.025624",
     "status": "completed"
    },
    "tags": []
   },
   "source": [
    "## 7. uMAP"
   ]
  },
  {
   "cell_type": "markdown",
   "id": "a216259f",
   "metadata": {
    "papermill": {
     "duration": 0.00351,
     "end_time": "2026-01-30T04:26:46.038153",
     "exception": false,
     "start_time": "2026-01-30T04:26:46.034643",
     "status": "completed"
    },
    "tags": []
   },
   "source": [
    "## 8. t-SNE"
   ]
  }
 ],
 "metadata": {
  "kaggle": {
   "accelerator": "none",
   "dataSources": [],
   "dockerImageVersionId": 30749,
   "isGpuEnabled": false,
   "isInternetEnabled": true,
   "language": "r",
   "sourceType": "notebook"
  },
  "kernelspec": {
   "display_name": "R",
   "language": "R",
   "name": "ir"
  },
  "language_info": {
   "codemirror_mode": "r",
   "file_extension": ".r",
   "mimetype": "text/x-r-source",
   "name": "R",
   "pygments_lexer": "r",
   "version": "4.4.0"
  },
  "papermill": {
   "default_parameters": {},
   "duration": 39.129877,
   "end_time": "2026-01-30T04:26:46.365752",
   "environment_variables": {},
   "exception": null,
   "input_path": "__notebook__.ipynb",
   "output_path": "__notebook__.ipynb",
   "parameters": {},
   "start_time": "2026-01-30T04:26:07.235875",
   "version": "2.6.0"
  }
 },
 "nbformat": 4,
 "nbformat_minor": 5
}
