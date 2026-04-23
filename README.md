# 🧬 Gene Analysis Toolkit

Analyze your raw DNA data with modular scripts for ancestry, disease risk, longevity, and fitness insights.

---

## 🚀 Quick Start

1. Place your genome file in the project folder  
2. Rename it to `Genome.txt`  
3. Run:

```bash
python run_all_analyses.py "/Users/yourname/Desktop/Genome.txt"
```

---

## 📦 Features

### 🌍 Ancestry / Ethnicity
- Population markers  
- Trait breakdown  
- Interpretable summaries  

### 🧠 Disease Risk
- SNP-based risk analysis  
- Polygenic risk scoring  
- Clear, simplified outputs  

### ⏳ Longevity
- Aging-related markers  
- Protective vs. risk variants  

Scripts:
- `Longevity.py` → quick summary  
- `Comprehensive_Longevity.py` → full probability analysis  

### 💪 Fitness
- Athletic trait probabilities  
- Performance-related genetic insights  

---

## 📄 Input Format

All scripts require a file named:

```
Genome.txt
```

### File Requirements
- Tab-separated text file  
- 4 columns (in order):

```
rsid    chromosome    position    genotype
```

- Lines starting with `#` are ignored (comments)

### Example

```
# rsid    chromosome  position  genotype
rs3094315	1	742429	AA
rs12562034	1	758311	GG
rs3934834	1	995669	CT
```

---

## 📥 Getting Your Genome File

If you’re using a personal genomics service:

1. Log in  
2. Go to **Settings → DNA Relatives**  
3. Click **Download Raw Data**  
4. Rename the file to `Genome.txt`  

---

## ⚠️ Disclaimer

- This tool is for **educational purposes only**  
- It **does not provide medical diagnoses**  
- Results are based on currently known genetic associations  
- Data is biased toward populations of European ancestry and may not generalize  

---

## 🧠 Notes

- Not all genetic variants are known or included  
- Environmental and lifestyle factors are **not accounted for**  
- Interpret results cautiously and consult professionals for medical decisions  
