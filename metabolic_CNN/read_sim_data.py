import json
from pathlib import Path
import pandas as pd
import numpy as np
import random
from sklearn.preprocessing import StandardScaler

class read_sim_data():
    def __init__(self, config_dir: Path) -> None:
        """
        Object to train and test metabolic CNN given a config json file.
        """
        # Make sure config_dir is a Path
        if not isinstance(config_dir, Path):
            config_dir = Path(config_dir)
        # Read config.json file
        with open(config_dir, 'r') as j:
            config_dict = json.loads(j.read())
            
        self.sim_data_dir = Path(config_dict['sim_data_dir']) / (config_dict['sim_data_filename'])
        self.sim_param_dir = Path(config_dict['sim_data_dir']) / (config_dict['sim_param_filename'])
            
        self.max_time = config_dict['max_time'] if not (config_dict['max_time'] == '') else 4
        self.min_time = config_dict['min_time'] if not (config_dict['min_time'] == '') else 2
        self.delta_t = config_dict['delta_t_sim'] if not (config_dict['delta_t_sim'] == '') else 0.1
        
        self.cols_sim_data = config_dict['cols_sim_data']
        self.cols_sim_param = config_dict['cols_sim_param']
        self.order_in_mids = config_dict['order_input_mids']
        
        self.n_metabs = config_dict['n_metabs']
        self.n_mids = config_dict['n_mids']
        self.n_t = int((self.max_time - self.min_time) / self.delta_t + 1)
        
        self.metab_sources = config_dict['metab_sources']
        self.target_flux = config_dict['target_flux']
        self.name_rel_flux = config_dict['name_rel_flux']
        self.non_mid_cols = config_dict['non_mid_cols']
        
        self.combined_data = self._read_sim_data()
        self.rel_target_flux = self._calc_rel_target_flux()
        
        self.combined_data = self.combined_data.drop(columns=self.metab_sources, 
                                                     axis=1)
        
        self.hypo_mids = config_dict['hypothetical_MIDs']
        if self.hypo_mids != '':
            for i in range(len(self.hypo_mids)):
                self.combined_data[self.hypo_mids[i]] = 0
                
        self.combined_data[self.name_rel_flux] = self.rel_target_flux
        
        config_dict['non_mid_cols'].append(self.name_rel_flux)
        self.non_mid_data = config_dict['non_mid_cols']
        self._rm_neg_mid()
        
        
    def _read_sim_data(self) -> pd.DataFrame:
        """
        Reads sim_data_filename and sim_param_filename defined in config file.
        Simulated data include index of simulation, time points, and balanced 
        metabolite MIDs. 
        Simulated parameters include input metabolite MIDs and fluxes of 
        sources of a metabolite.

        Returns an integrated dataframe of simulated data and simulated 
        parameters.
        Rows of combined_mid_param: number of simulations * number of time points
        columns of combined_mid_param: index, time, MIDs of related metabolites,
        source fluxes
        """
        sim_data = pd.read_csv(self.sim_data_dir)
        sim_param = pd.read_csv(self.sim_param_dir)
        
        sim_data = sim_data.loc[sim_data['time'] >= self.min_time] 
        sim_data = sim_data[self.cols_sim_data]
        sim_param = sim_param[self.cols_sim_param]        
        
        self.n_sim = len(sim_param)
        
        df_temp_all = []
        for i in range(self.n_sim):
            df_temp = pd.concat([sim_param.iloc[i]] * self.n_t, axis = 1)
            df_temp = df_temp.T
            df_temp_all.append(df_temp)

        df_temp_all = pd.concat(df_temp_all)

        combined_mid_param = pd.concat([sim_data.reset_index(drop=True), 
                                        df_temp_all.reset_index(drop=True)], 
                                       axis = 1)
        ord_params = self.non_mid_cols + self.order_in_mids + self.metab_sources
        combined_mid_param = combined_mid_param.reindex(ord_params, axis = 1)
        
        return combined_mid_param
    
    def _calc_rel_target_flux(self) -> pd.Series:
        """
        Calculates ratio of target flux over all sources of a metabolite

        Returns a pd.Series of a relative target flux with the shape of number of
        simulations * number of time points
        """
        sum_source_flux = np.zeros(self.combined_data.shape[0])
        for i in range(len(self.metab_sources)):
            sum_source_flux = sum_source_flux + self.combined_data[self.metab_sources[i]]
            
        rel_target_flux = self.combined_data[self.target_flux] / sum_source_flux

        return rel_target_flux
    
    def _reshape_nsim_mid(self, data_df: pd.DataFrame, n_sim: int) -> np.ndarray:
        """
        Reshapes feature matrix (i.e., MIDs) to (number of simulations, 
                                                 number of time points * 
                                                 number of metabolites * 
                                                 number of isotopologues)

        Args:
            data_df (pd.DataFrame): a feature matrix including balanced and input MIDs
        
            n_sim (int): number of simulation instances
        
        Returns:
            data_np (np.ndarray): reshaped feature matrix with n_sim rows
        """
        data_np = data_df.drop(columns = self.non_mid_data, axis = 1)             
        data_np = data_np.to_numpy()
        data_np = np.reshape(data_np, 
                             (n_sim, self.n_t * self.n_metabs * self.n_mids))
        return data_np 
    
    def _rm_neg_mid(self) -> None:
        """
        Removes simulations that generate negative MIDs if any.
        """
        nsim_mid = self._reshape_nsim_mid(self.combined_data, self.n_sim)
        
        neg_flux = np.sum(nsim_mid < 0, axis = 1)
        neg_flux_ind = np.nonzero(neg_flux)
        nsim_mid = np.delete(nsim_mid, neg_flux_ind, axis = 0)    
        
        self.n_sim = nsim_mid.shape[0]
        
        self.combined_data = self.combined_data.loc[~self.combined_data['index'].isin(neg_flux_ind[0])]
        
    def split_sim_data(self, data_df: pd.DataFrame, split_frac = 0.15) -> tuple:
        """
        Shuffles and splits data into two subsets. 

        Args:
            data_df (pd.DataFrame): an integrated dataframe of simulated data 
            and simulated parameters.
            
            split_frac (float): a fraction that splits data into the first 
            dataset (1-test_frac) and the second dataset (test frac).The 
            default is 0.15.

        Returns a tuple of (first subset, second subset)
        """
        # shuffle across simulations
        sim_index = data_df['index'].unique() 
        random.seed(123)
        sub2_index = random.sample(list(sim_index), round(len(sim_index) * split_frac))
    
        sub2_data = data_df.loc[data_df['index'].isin(sub2_index)]
        sub1_data = data_df.loc[~(data_df['index'].isin(sub2_index))] 
        
        return(sub1_data, sub2_data)
    
    def fit_scaler(self, train_data: pd.DataFrame):
        """
        Computes mean and standard deviation of MIDs across training simulations

        Args:
            train_data (pd.DataFrame): an integrated dataframe of simulated data 
            and simulated parameters for training 

        Returns:
            scaler (object): fitted scaler
        """
        nsim_train = len(np.unique(train_data['index'].to_numpy()))
        train_data = self._reshape_nsim_mid(train_data, nsim_train)
        
        scaler = StandardScaler()
        scaler.fit(train_data)
        return scaler
    
    def transform_scaler_sim_data(self, scaler, data_df: pd.DataFrame) -> pd.DataFrame:
        """
        Standardize MIDs by centering and scaling

        Args:
            scaler (object): fitted scaler
            data_df (pd.DataFrame): an integrated dataframe of simulated data 
            and simulated parameters

        Returns:
            data_df (pd.DataFrame): transformed dataframe
        """
        nsim_data = len(np.unique(data_df['index'].to_numpy()))
        data_np = self._reshape_nsim_mid(data_df, nsim_data)
                
        data_np = scaler.transform(data_np)
        data_np = np.reshape(data_np, (nsim_data * self.n_t, 
                                       self.n_metabs * self.n_mids))
        
        #replace transformed mids into original dataframe
        data_df.iloc[:, 2:2 + self.n_metabs * self.n_mids] = data_np
        
        return data_df
    
    def split_features_target(self, data_df: pd.DataFrame) -> tuple:
        """
        Splits features (i.e., MIDs) from target (i.e., relative target flux)
        
        Args:
            data_df (pd.DataFrame): an integrated dataframe of simulated data 
            and simulated parameters 
            
        Returns a tuple of (features, target)
        """
        # shuffle the DataFrame rows
        data_df = data_df.sample(frac = 1, random_state = 123)
    
        # Divide the data into features and target
        features = data_df.drop(columns=self.non_mid_data)
        target = data_df[self.name_rel_flux].to_numpy()
        features = features.to_numpy()
        
        return(features, target)
    
    def reshape_features(self, data_np: np.ndarray) -> np.ndarray:
        """
        Aligns features (i.e., MIDs) to the shape of CNN input layer

        Args:
            data_np (np.ndarray): features with the shape of (number of simulations *
                                                              number of time points,
                                                              number of metabolites *
                                                              number of isotopologues)

        Returns:
            data_np (np.ndarray): features with the shape of (number of simulations *
                                                              number of timepoints,
                                                              number of metabolites,
                                                              number of isotopologues)
        """
        data_np = np.reshape(data_np, (data_np.shape[0], 
                                       self.n_metabs, self.n_mids, 1))
        return data_np
    