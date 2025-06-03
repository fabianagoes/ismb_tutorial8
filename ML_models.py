from sklearn.datasets import load_iris
from sklearn.model_selection import train_test_split
from sklearn.metrics import accuracy_score
from sklearn.tree import DecisionTreeClassifier
from sklearn.ensemble import RandomForestClassifier
from sklearn.linear_model import LogisticRegression
from sklearn.neural_network import MLPClassifier
from xgboost import XGBClassifier
from tqdm import tqdm

models = {
    "Decision Tree": DecisionTreeClassifier(),
    "Random Forest": RandomForestClassifier(),
    "Logistic Regression": LogisticRegression(max_iter=1000),
    #"MLP (Neural Net)": MLPClassifier(max_iter=1000),
    "XGBoost": XGBClassifier(use_label_encoder=False, eval_metric='mlogloss')
}

def run_models(Xtrain, Ytrain, Xtest, Ytest):
    acc=[]
    for name in tqdm(models, desc="Training Models"):
        clf = models[name]
        clf.fit(Xtrain,Ytrain)
        y_pred = clf.predict(Xtest)
        acc.append(accuracy_score(Ytest, y_pred))
    return acc