import numpy as np
import matplotlib.pyplot as plt
import src.utils_energy_distribution as en_util


def display_optimization_stats(loss, loss_per_epoch, lr, fig_width=10):
    # Data for plotting
    iters_range = np.arange(len(loss)) + 1
    epochs_range = np.arange(len(loss_per_epoch)) + 1

    # Create figure
    fig, (ax1, ax2, ax3) = plt.subplots(3, 1)
    fig.set_figwidth(fig_width)
    fig.set_figheight(3 * fig_width)

    # loss function axis
    ax1.semilogy(iters_range, loss)
    ax1.set(title='Loss per iters')
    ax1.grid()

    # LR axis
    ax2.plot(iters_range, lr)
    ax2.set(title='LR per iters')
    ax2.grid()

    # Epoch loss axis
    ax3.plot(epochs_range, loss_per_epoch)
    ax3.set(title='Epoch Loss')
    ax3.grid()

    fig.tight_layout()
    plt.show()


def display_random_pairs(ref_funcs_list, meas_funcs_list, titles=[], i=-1):
    num_pairs = len(ref_funcs_list)
    num_funcs = ref_funcs_list[0].shape[0]
    if i < 0:
        i = np.random.randint(0, num_funcs)

    if len(titles) < num_pairs * 2:
        do_defaul_titles = True
    else:
        do_defaul_titles = False

    display_list = []

    for j in range(num_pairs):
        ref_func = ref_funcs_list[j][i]
        meas_func = meas_funcs_list[j][i]
        rrmse = en_util.masked_rrmse_loss(ref_func, meas_func)

        display_list.append(ref_func)
        display_list.append(meas_func)

        if do_defaul_titles:
            titles.append(str(j) + ' - ref')
            titles.append(str(j) + ' - meas')

        k = j * 2 + 1
        titles[k] = titles[k] + f" (RRMSE={rrmse :.2f}%)"

    print(f"Index={i :d}")
    display_figs(display_list, titles)
    return i


def display_2D_func(func_vals, arg_vals=[], fig_width=10, title='2D Plot'):
    if not arg_vals:
        arg_vals = np.arange(len(func_vals)) + 1
    fig, ax = plt.subplots(1, 1)
    fig.set_figwidth(fig_width)
    fig.set_figheight(fig_width)
    ax.plot(arg_vals, func_vals)
    ax.set(title=title)
    ax.grid()
    plt.show()


def display_figs(display_list, titles, fig_height=7):
    plt.figure(figsize=(fig_height * len(display_list), fig_height))

    for i in range(len(display_list)):
        plt.subplot(1, len(display_list), i + 1)
        plt.title(titles[i])
        plt.imshow(display_list[i])
        plt.axis('off')
    plt.show()
