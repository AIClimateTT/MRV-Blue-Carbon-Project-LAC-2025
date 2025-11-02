# MRV-Blue-Carbon-Project-LAC-2025 - CARONI PILOT
The goal of the Caroni pilot is to produce **accurate, scalable, and reproducible estimates of mangrove above-ground carbon (AGC)** for the **Caroni Swamp in Trinidad and Tobago**, as part of the broader regional AI MRV Blue Carbon initiative. This work directly supports the country’s Measurement, Reporting, and Verification (MRV) obligations under international climate agreements by developing a locally calibrated, AI-driven model that can be replicated and scaled across the Caribbean.

The model produces pixel-level AGC estimates (Mg per 20 m) for the Caroni mangrove extent and validates results against 2018 field plot data and LiDAR-derived carbon surfaces. These estimates serve as a technical proof of concept for scaling the workflow to other countries (Colombia, Jamaica, Panama, Suriname) under the AI MRV Blue Carbon framework.

Although this phase focuses primarily on above-ground carbon, the resulting workflow is designed to integrate Soil Organic Carbon (SOC) and methane flux data as these become available. This ensures future compatibility with results-based carbon finance mechanisms and carbon credit market verification standards.

# Caroni Model Architecture & Data Integration

The Caroni pipeline applies a supervised machine learning approach (Random Forest, XGBoost) using spectral, structural, and environmental predictors derived from satellite data and auxiliary sources:

1) Predictors: Sentinel-2 bands (B02–B12), NDVI, and optional DEM values.
2) Target variable: derived AGC (Mg per 20 m pixel).
3) Training data: UWI-provided AGC tables (2015–2025), validated against 2018 field plots.
4) Spatial extent: Caroni Swamp mangrove area (UTM Zone 20N, ~108,000 pixels).
5) Output products: Continuous AGC surface raster (20 m), feature importance analysis, and validation statistics (R², MAE, MSE).

All models are developed and tested within a reproducible Colab/Python environment, using version-controlled scripts and standardized metadata conventions.

Repository Structure (Caroni Pilot)

All scripts, data preparation code, and metadata for the Caroni pilot are organized in the AI MRV Blue Carbon – Caroni Repository, following open-source, reproducible data science standards.

Directory	Description
/01_raw/	Original input data (Sentinel-2 Excel files, LiDAR AGC raster, field plots).
/02_interim/	Cleaned and co-registered tables; merged pixel-level datasets (AGC + bands + NDVI + DEM).
/03_features/	Derived features and vegetation indices (NDVI, red-edge indices, etc.).
/04_models/	Model training notebooks, parameter logs, serialized RF/XGB models.
/05_outputs/	AGC prediction rasters, uncertainty layers (P10/P50/P90), validation metrics, and feature-importance charts.

# Regional Transfer Learning

The Caroni pilot serves as the baseline model for regional transfer learning within the AI MRV Blue Carbon initiative. By establishing a fully validated end-to-end workflow—covering data preprocessing, feature extraction, and machine learning model calibration. The Caroni model provides a reference framework that can be fine-tuned using country-specific data in Colombia, Jamaica, Panama, and Suriname. This approach leverages shared spectral–ecological patterns across mangrove ecosystems while allowing for local calibration through additional field plots or LiDAR acquisitions. The Caroni model thus acts as a foundational “source model” for regional upscaling, reducing computational and data preparation costs, and ensuring methodological consistency across all participating countries.
