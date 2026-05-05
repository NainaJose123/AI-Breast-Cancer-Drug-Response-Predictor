🧬 AI-Based Breast Cancer Drug Response Predictor

A bioinformatics-driven machine learning system that predicts breast cancer patient response to targeted therapies using multi-omics gene expression data and pharmacogenomic datasets. The project integrates biological knowledge with AI to support precision medicine applications.




🚀 Overview
This system predicts how a patient’s tumor may respond to targeted breast cancer drugs such as:
Trastuzumab
Gefitinib
Lapatinib
It classifies drug response into:
🟥 Resistant
🟧 Moderate
🟩 Sensitive
The project combines machine learning, gene expression analysis, and biological pathway knowledge to enable personalized cancer therapy prediction.




🧠 Key Features
End-to-end ML pipeline for drug response prediction
Multi-drug support with configurable selection
High-dimensional gene expression processing
Real-time prediction via FastAPI backend
CSV-based patient input system
Interactive web-based simulation interface
Extensible cascading decision architecture (under refinement)
Integration of biological pathway knowledge for interpretability




🧬 Data Sources & Biological Databases


This project integrates multiple real-world biomedical datasets and knowledge bases:


Gene Expression Data:
GEO GSE81538
Used for tumor RNA-seq gene expression profiles.


Drug Response Data:
GDSC v2
Provides IC₅₀ drug sensitivity measurements.


Cancer Cell Line Data:
CCLE
Links gene expression with pharmacological response.


Protein Interaction Network:
STRING v12.0
Used for modeling functional gene/protein interactions.


Drug–Target Knowledge Base:
DrugBank
Used to map drugs to biological targets (e.g., ERBB2, PIK3CA).


Pathway Information:
KEGG pathways
Used for cancer pathway enrichment and biological validation.




⚙️ Methodology
1. Data Preprocessing
RNA-seq gene expression cleaning
Normalization and scaling
Missing value handling

2. Feature Engineering
High-dimensional gene filtering
Selection of biologically relevant features
Integration with drug-target interaction data

3. Machine Learning Model
Algorithm: Random Forest Classifier
Output: Drug response category
Designed for extensibility to other ML models

4. Biological Integration
Drug-target mapping using DrugBank
Pathway alignment using KEGG
Interaction context using STRING network




🧩 System Architecture
FastAPI backend for prediction service
Modular ML pipeline:
Data preprocessing → Feature engineering → Prediction
CSV-based batch inference support
API-driven frontend interface
Designed for cascading multi-stage decision flow (in progress)




🌐 Web Interface
The system provides an interactive interface supporting:
Upload of patient gene expression CSV files
Drug selection for prediction
Tumor subtype filtering
Real-time prediction results via API




📊 API Endpoints
GET /ui → Web-based interface
POST /predict → Drug response prediction API




🔧 Installation & Setup
pip install -r requirements.txt
Run backend server:
uvicorn main:app --reload
Open in browser:
http://127.0.0.1:8000/ui




⚠️ Current Limitations
Cascading multi-stage decision pipeline is under refinement
Stage-wise data flow consistency needs improvement
Model optimization for large-scale genomic data is ongoing




🔮 Future Improvements
Integration of deep learning models (CNN/Transformer for omics data)
Improved cascading decision system with explainability
SHAP-based feature importance analysis
Cloud deployment (Render / AWS / HuggingFace Spaces)
Addition of survival prediction module
Enhanced visualization of gene–drug interaction networks




🧑‍💻 Tech Stack
Python
FastAPI
Scikit-learn
Pandas / NumPy
Bioinformatics datasets (GEO, GDSC, CCLE, STRING, DrugBank, KEGG)




📌 Impact
This project bridges machine learning and cancer biology by integrating multi-source biomedical datasets to support precision medicine and drug response prediction in breast cancer.




🖼️ System Preview
<img width="719" height="305" alt="image" src="https://github.com/user-attachments/assets/3f1894c9-41ec-4bd0-b42b-ca9120d9e5b2" />
<img width="715" height="255" alt="image" src="https://github.com/user-attachments/assets/80a08da5-13b9-49aa-a521-ca3a88aec00d" />
<img width="719" height="318" alt="image" src="https://github.com/user-attachments/assets/160e2a91-c512-4c20-98a8-9159a901fffb" />
<img width="719" height="283" alt="image" src="https://github.com/user-attachments/assets/3e2e8824-4b96-47c0-a11a-2c1ec9823c11" />
<img width="714" height="275" alt="image" src="https://github.com/user-attachments/assets/c873875f-d804-42a7-993a-2b2073a51d52" />
<img width="714" height="341" alt="image" src="https://github.com/user-attachments/assets/13dd262c-8f5e-47a9-843f-de12b15e388a" />








🎥 Project Demo
https://drive.google.com/file/d/1SON7ZbocyueQ-RfFAAJvCIyAh2fHuXM6/view?usp=drive_link




👤 Author
Developed as an AI + Bioinformatics research project focused on personalized cancer therapy prediction using multi-omics data integration.
