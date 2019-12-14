import pandas
import numpy as np
import sklearn
from sklearn.svm import SVC
from sklearn.model_selection import KFold, GridSearchCV


def select_table(table, sets):
    """ Select a subset of table with certain genes
    Args:
        table: pandas dataframe with rows are samples, columns are genes
        sets: list of sets in which each set contains interested genes
    Returns:
        subset table
    """
    
    union_set = set()
    for s in sets:
        union_set = union_set.union(s)
    return table[list(union_set)]
    


def training(X_train, y_train):
    """ Perform cross_validation on SVM
    Args:
        X_train: input data, expected to has numpy format with rows are samples, columns are genes
        y_train: output
    Returns:
        Trained model
    """

    possible_parameters = {
        'C': [1e1, 1e2, 1e3],
        'gamma': [1e0, 1e-1, 1e-2]
    }
    svc = SVC(kernel='rbf')
    
    # Perform grid search cross validation
    clf = GridSearchCV(svc, possible_parameters, n_jobs=8, cv=3)
    clf.fit(X_train, y_train)
    print(clf.best_params_)
    
    # Train with best setting to get scores on cross validation
    kf = KFold(n_splits=3, shuffle=True)
    scores = sklearn.model_selection.cross_validate(clf.best_estimator_, X_train, y_train, cv=kf, scoring=['accuracy', 'f1', 'precision', 'recall'])
    for key in scores:
        print(key, np.mean(scores[key]))
        
    return clf.best_estimator_
