from pathlib import Path

import ase
import matplotlib.pyplot as plt
import numpy as np
import optax
from dcmnet.analysis import create_model_and_params
from dcmnet.utils import apply_model
from tqdm import tqdm


def atom_centered_dipole(dcm, com, q):
    dipole_out = np.zeros(3)
    for i, _ in enumerate(dcm):
        dipole_out += q[i] * (_ - com)
    # print(dipole_out*2.5417464519)
    return np.linalg.norm(dipole_out) * 4.80320


batch_size = 1


ddir = Path("/pchem-data/meuwly/boittier/home/jaxeq/all_runs/")
files = {_.parents[0].name: _ for _ in list(ddir.glob("dcr175/*dcm*/best*"))}
KEYS = list(files.keys())


def get_dcmnet_predictions(batch, verbose=False):
    q_predictions = []
    qr_predictions = []
    rmse_predictions = []
    dipole_predictions = []
    N_dcms = []

    for i in range(len(KEYS)):
        p = files[KEYS[i]]
        model, params, _ = create_model_and_params(p, debug=False)
        if verbose:
            print(model)
            print(p)
        mask = batch["espMask"]
        mono, dipo = apply_model(model, params, batch, 1)

        print(mono.shape, dipo.shape)

        from dcmnet.loss import esp_mono_loss_pots

        if model.n_dcm == 1:
            dipo = batch["R"]

        esp_dc_pred = esp_mono_loss_pots(
            dipo,
            mono.flatten(),
            batch["vdw_surface"],
            batch["Z"],
            batch_size,
            model.n_dcm,
        )

        plt.scatter(batch["esp"][mask], esp_dc_pred[0][mask])
        ax = plt.gca()
        ax.set_xlim(-0.1, 0.1)
        ax.set_ylim(-0.1, 0.1)
        ax.set_aspect("equal")
        l2_loss = optax.l2_loss(batch["esp"][mask], esp_dc_pred[0][mask]) * 2
        esp_loss = np.mean(l2_loss) ** 0.5
        pred_D = atom_centered_dipole(dipo, batch["com"], mono.flatten())
        dipole_str = "DIPOLE (DEBYE): {:.3f}".format(pred_D)
        rmse_str = "RMSE ESP: {:.3f}".format(esp_loss * 627.509)
        plt.title(f"{dipole_str}\n{rmse_str}\n N DCM: {model.n_dcm}")
        plt.show()

        q_predictions.append(mono)
        qr_predictions.append(dipo)
        dipole_predictions.append(pred_D)
        rmse_predictions.append(esp_loss * 627.509)
        N_dcms.append(mono.shape[2] * batch["N"][0])

        dipo = dipo.reshape(mono.shape[1], mono.shape[2], 3)

        # now sum the hydrogen charges and average their positions
        new_mono = []
        new_dipo = []
        for i in range(len(mono[0])):
            if batch["Z"][i] == 1:
                print(i, "hydrogen sum", mono[0][i].shape, dipo[i].shape)
                new_mono.append(mono[0][i].sum())
                new_dipo.append(dipo[i].mean(axis=-2))
            elif batch["Z"][i] != 0:
                for j in range(len(mono[0][i])):
                    new_mono.append(mono[0][i][j])
                    new_dipo.append(dipo[i][j])

        mono = np.array(new_mono)
        dipo = np.array(new_dipo).reshape(len(mono), 3)
        print("hydrogen sum", mono.shape, dipo.shape)
        print("hydrogen sum", mono, dipo)

        q_predictions.append(mono)
        qr_predictions.append(dipo)
        dipole_predictions.append(pred_D)
        rmse_predictions.append(esp_loss * 627.509)
        N_dcms.append(mono.shape[0])

    # q_predictions, qr_predictions, rmse_predictions, dipole_predictions, N_dcms
    out_dict = {
        "q_predictions": q_predictions,
        "qr_predictions": qr_predictions,
        "rmse_predictions": rmse_predictions,
        "dipole_predictions": dipole_predictions,
        "N_dcms": N_dcms,
    }
    return out_dict
